package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import htsjdk.samtools.SAMFlag;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFastqUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndexCache;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembler;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembly;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import scala.Tuple2;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.*;

/** This LocalAssemblyHandler aligns assembly contigs with BWA, along with some optional writing of intermediate results. */
public final class FermiLiteAssemblyHandler implements FindBreakpointEvidenceSpark.LocalAssemblyHandler {
    private static final long serialVersionUID = 1L;
    private final String alignerIndexFile;
    private final int maxFastqSize;
    private final String fastqDir;
    private final boolean writeGFAs;

    public FermiLiteAssemblyHandler( final String alignerIndexFile, final int maxFastqSize,
                                     final String fastqDir, final boolean writeGFAs ) {
        this.alignerIndexFile = alignerIndexFile;
        this.maxFastqSize = maxFastqSize;
        this.fastqDir = fastqDir;
        this.writeGFAs = writeGFAs;
    }

    @Override
    public AlignedAssemblyOrExcuse apply( final Tuple2<Integer, List<SVFastqUtils.FastqRead>> intervalAndReads ) {
        final int intervalID = intervalAndReads._1();
        final String assemblyName = AlignedAssemblyOrExcuse.formatAssemblyID(intervalID);
        final List<SVFastqUtils.FastqRead> readsList = intervalAndReads._2();

        final int fastqSize = readsList.stream().mapToInt(FastqRead -> FastqRead.getBases().length).sum();
        if ( fastqSize > maxFastqSize ) {
            return new AlignedAssemblyOrExcuse(intervalID, "no assembly -- too big (" + fastqSize + " bytes).");
        }

        if ( fastqDir != null ) {
            final String fastqName = String.format("%s/%s.fastq", fastqDir, assemblyName);
            final ArrayList<SVFastqUtils.FastqRead> sortedReads = new ArrayList<>(readsList);
            sortedReads.sort(Comparator.comparing(SVFastqUtils.FastqRead::getHeader));
            SVFastqUtils.writeFastqFile(fastqName, sortedReads.iterator());
        }

        final long timeStart = System.currentTimeMillis();
        final FermiLiteAssembly assembly = new FermiLiteAssembler().createAssembly(readsList);
        final int secondsInAssembly = (int)((System.currentTimeMillis() - timeStart + 500)/1000);

        if ( fastqDir != null && writeGFAs ) {
            final String gfaName =  String.format("%s/%s.gfa", fastqDir, assemblyName);
            try ( final OutputStream os = BucketUtils.createFile(gfaName) ) {
                assembly.writeGFA(os);
            }
            catch ( final IOException ioe ) {
                throw new GATKException("Can't write "+gfaName, ioe);
            }
        }

        final List<byte[]> tigSeqs =
                assembly.getContigs().stream()
                        .map(FermiLiteAssembly.Contig::getSequence)
                        .collect(SVUtils.arrayListCollector(assembly.getNContigs()));
        final AlignedAssemblyOrExcuse result;
        try ( final BwaMemAligner aligner = new BwaMemAligner(BwaMemIndexCache.getInstance(alignerIndexFile)) ) {
            aligner.setIntraCtgOptions();
            final List<List<BwaMemAlignment>> alignments = aligner.alignSeqs(tigSeqs);
            result = new AlignedAssemblyOrExcuse(intervalID, assembly, secondsInAssembly, alignments);
        }
        final int nTigs = tigSeqs.size();
        if ( nTigs == 0 ) {
            return result;
        }

        final String tmpDir = "/tmp/fbes";
        try {
            IOUtils.createDirectory(tmpDir);
        } catch ( final IOException ioe ) {
            throw new GATKException("Can't create " + tmpDir, ioe);
        }
        final String fastaFile = String.format("%s/%s.fasta", tmpDir, assemblyName);
        try ( final BufferedOutputStream os = new BufferedOutputStream(BucketUtils.createFile(fastaFile)) ) {
            for (int tigId = 0; tigId != nTigs; ++tigId) {
                final String seqName = assemblyName + "." + tigId;
                final byte[] line1 = (">" + seqName + "\n").getBytes();
                final byte[] seq = tigSeqs.get(tigId);
                os.write(line1);
                os.write(seq);
                os.write('\n');
            }
        }
        catch ( final IOException ioe ) {
            throw new GATKException("Unable to write fasta file of contigs from assembly " + assemblyName, ioe);
        }

        final String imageFile = String.format("%s/%s.img", tmpDir, assemblyName);
        BwaMemIndex.createIndexImageFromFastaFile(fastaFile, imageFile);
        try { BucketUtils.deleteFile(fastaFile); }
        catch ( final IOException ioe ) { throw new GATKException("unable to delete " + fastaFile, ioe); }

        final List<List<BwaMemAlignment>> alignments;
        try ( final BwaMemIndex assemblyIndex = new BwaMemIndex(imageFile);
              final BwaMemAligner aligner = new BwaMemAligner(assemblyIndex) ) {
            aligner.alignPairs();
            alignments =
                    aligner.alignSeqs(readsList, SVFastqUtils.FastqRead::getBases);
        }

        try { BucketUtils.deleteFile(imageFile); }
        catch ( final IOException ioe ) { throw new GATKException("unable to delete " + imageFile, ioe); }

        final Map<Link, Integer> linkCounts = new HashMap<>();
        final int nReads = readsList.size();
        for ( int idx = 0; idx < nReads; idx += 2 ) {
            final List<BwaMemAlignment> alignList1 = alignments.get(idx);
            final List<BwaMemAlignment> alignList2 = alignments.get(idx + 1);
            for ( final BwaMemAlignment alignment1 : alignList1 ) {
                for ( final BwaMemAlignment alignment2 : alignList2 ) {
                    boolean isRC1 = SAMFlag.READ_REVERSE_STRAND.isSet(alignment1.getSamFlag());
                    boolean isRC2 = !SAMFlag.READ_REVERSE_STRAND.isSet(alignment2.getSamFlag());
                    if ( alignment1.getRefId() != alignment2.getRefId() || isRC1 != isRC2 ) {
                        final Link link =
                                new Link(new ContigStrand(alignment1.getRefId(), isRC1),
                                         new ContigStrand(alignment2.getRefId(), isRC2));
                        linkCounts.merge(link, 1, Integer::sum);
                    }
                }
            }
        }
        final Iterator<Map.Entry<Link, Integer>> linkItr = linkCounts.entrySet().iterator();
        while ( linkItr.hasNext() ) {
            final Map.Entry<Link, Integer> entry = linkItr.next();
            if ( entry.getValue() < 10 ) linkItr.remove();
            final Link link = entry.getKey();
            final FermiLiteAssembly.Contig sourceContig = assembly.getContig(link.getId1());
            final FermiLiteAssembly.Contig targetContig = assembly.getContig(link.getId2());
            boolean foundConnection = false;
            for ( final FermiLiteAssembly.Connection connection : sourceContig.getConnections() ) {
                if ( connection.getTarget() == targetContig &&
                        connection.isRC() == link.isRC1() &&
                        connection.isTargetRC() == link.isRC2() ) {
                    entry.setValue(connection.getOverlapLen());
                    foundConnection = true;
                    break;
                }

            }
            if ( !foundConnection ) linkItr.remove();
        }

        return result;
    }

    private static final class ContigStrand {
        private final int contigId;
        private final boolean isRC;

        public ContigStrand( final int contigId, final boolean isRC ) {
            this.contigId = contigId;
            this.isRC = isRC;
        }

        public int getId() { return contigId; }
        public boolean isRC() { return isRC; }

        @Override public boolean equals( final Object obj ) {
            return obj instanceof ContigStrand && equals((ContigStrand) obj);
        }

        public boolean equals( final ContigStrand that ) {
            return contigId == that.contigId && isRC == that.isRC;
        }

        @Override public int hashCode() {
            int hashVal = 113;
            hashVal = 47 * (hashVal + contigId);
            return 47 * (hashVal + (isRC ? 3 : 5));
        }
    }

    private static final class Link {
        private final ContigStrand contigStrand1;
        private final ContigStrand contigStrand2;

        /** Canonicalizing constructor: id1 <= id2 */
        public Link( final ContigStrand contigStrand1, final ContigStrand contigStrand2 ) {
            if ( contigStrand1.getId() < contigStrand2.getId() ) {
                this.contigStrand1 = contigStrand1;
                this.contigStrand2 = contigStrand2;
            } else {
                this.contigStrand1 = contigStrand2;
                this.contigStrand2 = contigStrand1;
            }
        }

        public int getId1() { return contigStrand1.getId(); }
        public boolean isRC1() { return contigStrand1.isRC(); }
        public int getId2() { return contigStrand2.getId(); }
        public boolean isRC2() { return contigStrand2.isRC(); }

        @Override public boolean equals( final Object obj ) {
            return obj instanceof Link && equals((Link)obj);
        }

        public boolean equals( final Link link ) {
            return this == link ||
                    (contigStrand1.equals(link.contigStrand1) && contigStrand2.equals(link.contigStrand2));
        }

        @Override public int hashCode() {
            int hashVal = 113;
            hashVal = 47 * (hashVal + contigStrand1.hashCode());
            return 47 * (hashVal + contigStrand2.hashCode());
        }
    }

    private static final class Adjacencies {
        private final List<ContigStrand> predecessors = new LinkedList<>();
        private final List<ContigStrand> successors = new LinkedList<>();

        public boolean isSource( final boolean rc ) { return (rc ? successors: predecessors).isEmpty(); }

        public void addPredecessor( final ContigStrand contigStrand ) { predecessors.add(contigStrand); }
        public void addSuccessor( final ContigStrand contigStrand ) { successors.add(contigStrand); }
    }
}
