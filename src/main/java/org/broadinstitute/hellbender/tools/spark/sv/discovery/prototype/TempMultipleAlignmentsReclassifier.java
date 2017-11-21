package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.RDDUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.util.*;
import java.util.stream.Collectors;

// TODO: 10/30/17 temporary class later to be merged to AssemblyContigAlignmentSignatureClassifier.java after #3752 is in
final class TempMultipleAlignmentsReclassifier {

    @DefaultSerializer(AssemblyContigWithFineTunedAlignments.Serializer.class)
    static final class AssemblyContigWithFineTunedAlignments {
        private static final AlignedContig.Serializer contigSerializer = new AlignedContig.Serializer();
        final AlignedContig updatedContig;
        final List<String> insertionMappings;

        AssemblyContigWithFineTunedAlignments(final AlignedContig updatedContig,
                                              final List<String> insertionMappings) {
            this.updatedContig = updatedContig;
            this.insertionMappings = Utils.nonNull(insertionMappings);
        }

        AssemblyContigWithFineTunedAlignments(final Kryo kryo, final Input input) {
            updatedContig = contigSerializer.read(kryo, input, AlignedContig.class);
            final int insertionMappingsSize = input.readInt();
            insertionMappings = new ArrayList<>(insertionMappingsSize);
            if (insertionMappingsSize != 0) {
                for (int i = 0; i < insertionMappingsSize; ++i) {
                    insertionMappings.add(input.readString());
                }
            }
        }

        // after fine tuning, a contig may have no good alignment left, or only 1
        final boolean isInformative() {
            return updatedContig.alignmentIntervals.size() > 1;
        }

        final boolean hasIncomePicture() {
            Utils.validateArg(updatedContig.alignmentIntervals.size() > 2,
                    "assumption that input contig has more than 2 alignments is violated.\n" +
                            updatedContig.toString());

            final AlignmentInterval head = updatedContig.alignmentIntervals.get(0),
                                    tail = updatedContig.alignmentIntervals.get(updatedContig.alignmentIntervals.size()-1);

            // tail not resuming the "direction of flow" started by head
            if ( !head.referenceSpan.getContig().equals(tail.referenceSpan.getContig()) || head.forwardStrand != tail.forwardStrand)
                return true;

            // reference span anomaly
            final SimpleInterval referenceSpanHead = head.referenceSpan,
                                 referenceSpanTail = tail.referenceSpan;
            if (referenceSpanHead.contains(referenceSpanTail) || referenceSpanTail.contains(referenceSpanHead))
                return true;

            final boolean refSpanAnomaly =
                    updatedContig.alignmentIntervals
                            .subList(1, updatedContig.alignmentIntervals.size() - 1)
                            .stream().map(ai -> ai.referenceSpan)
                            .anyMatch( middle -> { // middle alignments' ref span either disjoint from head/tail, or completely contained in head/tail ref span
                                final boolean badHead = middle.overlaps(referenceSpanHead) && !referenceSpanHead.contains(middle);
                                final boolean badTail = middle.overlaps(referenceSpanTail) && !referenceSpanTail.contains(middle);
                                return badHead || badTail;
                            });
            if (refSpanAnomaly)
                return true;

            if (head.forwardStrand) {
                return referenceSpanHead.getStart() >= referenceSpanTail.getStart();
            } else {
                return referenceSpanHead.getEnd() <= referenceSpanHead.getEnd();
            }
        }

        void serialize(final Kryo kryo, final Output output) {
            contigSerializer.write(kryo, output, updatedContig);
            output.writeInt(insertionMappings.size());
            for (final String mapping : insertionMappings) {
                output.writeString(mapping);
            }
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<AssemblyContigWithFineTunedAlignments> {
            @Override
            public void write(final Kryo kryo, final Output output, final AssemblyContigWithFineTunedAlignments alignedContig) {
                alignedContig.serialize(kryo, output);
            }

            @Override
            public AssemblyContigWithFineTunedAlignments read(final Kryo kryo, final Input input, final Class<AssemblyContigWithFineTunedAlignments> clazz) {
                return new AssemblyContigWithFineTunedAlignments(kryo, input);
            }
        }
    }

    /**
     * convenience struct holding 4 classes:
     *      1) non-informative contigs who after fine tuning has 0 or 1 good alignment left
     *      2) contigs with more than 2 good alignments but doesn't seem to have picture complete as defined by
     *          {@link AssemblyContigWithFineTunedAlignments#hasIncomePicture()}
     *      3) contigs with 2 good alignments and bad alignments encoded as strings
     *      4) contigs with more than 2 good alignments and seemingly have picture complete as defined by
     *          {@link AssemblyContigWithFineTunedAlignments#hasIncomePicture()}
     */
    private static final class MultipleAlignmentReclassificationResults {
        final JavaRDD<AssemblyContigWithFineTunedAlignments> nonInformativeContigs;
        final JavaRDD<AssemblyContigWithFineTunedAlignments> contigsWithMoreThanTwoGoodAlignments;
        final JavaRDD<AssemblyContigWithFineTunedAlignments> contigsWithTwoGoodAlignments;
        final JavaRDD<AssemblyContigWithFineTunedAlignments> contigsWithMoreThanTwoGoodAlignmentsButIncompletePicture;
        MultipleAlignmentReclassificationResults(final JavaRDD<AssemblyContigWithFineTunedAlignments> nonInformativeContigs,
                                                 final JavaRDD<AssemblyContigWithFineTunedAlignments> contigsWithMoreThanTwoGoodAlignments,
                                                 final JavaRDD<AssemblyContigWithFineTunedAlignments> contigsWithTwoGoodAlignments,
                                                 final JavaRDD<AssemblyContigWithFineTunedAlignments> contigsWithMoreThanTwoGoodAlignmentsButIncompletePicture) {
            this.nonInformativeContigs = nonInformativeContigs;
            this.contigsWithMoreThanTwoGoodAlignments = contigsWithMoreThanTwoGoodAlignments;
            this.contigsWithTwoGoodAlignments = contigsWithTwoGoodAlignments;
            this.contigsWithMoreThanTwoGoodAlignmentsButIncompletePicture = contigsWithMoreThanTwoGoodAlignmentsButIncompletePicture;
        }
    }

    // =================================================================================================================

    /**
     * Reclassify assembly contigs based on alignment fine tuning as implemented in
     * {@link #removeNonUniqueMappings(List, int, int)}.
     */
    static MultipleAlignmentReclassificationResults reClassifyContigsWithMultipleAlignments(
            final JavaRDD<AlignedContig> localAssemblyContigs,
            final int mapQThresholdInclusive, final int uniqReadLenInclusive) {

        final JavaRDD<AssemblyContigWithFineTunedAlignments> contigsWithFineTunedAlignments =
                localAssemblyContigs
                        .map(tig -> {
                            final Tuple2<List<AlignmentInterval>, List<AlignmentInterval>> goodAndBadAlignments =
                                    removeNonUniqueMappings(tig.alignmentIntervals, mapQThresholdInclusive, uniqReadLenInclusive);
                            final List<AlignmentInterval> goodAlignments = goodAndBadAlignments._1;
                            final List<String> insertionMappings =
                                    goodAndBadAlignments._2.stream().map(AlignmentInterval::toPackedString).collect(Collectors.toList());

                            return
                                    new AssemblyContigWithFineTunedAlignments(
                                            new AlignedContig(tig.contigName, tig.contigSequence, goodAlignments,
                                                    tig.hasEquallyGoodAlnConfigurations),
                                            insertionMappings);
                        });

        // first take down non-informative assembly contigs
        final Tuple2<JavaRDD<AssemblyContigWithFineTunedAlignments>, JavaRDD<AssemblyContigWithFineTunedAlignments>> informativeAndNotSo =
                RDDUtils.split(contigsWithFineTunedAlignments, TempMultipleAlignmentsReclassifier::isInformative, false);
        final JavaRDD<AssemblyContigWithFineTunedAlignments> garbage = informativeAndNotSo._2;

        // assembly contigs with 2 good alignments and bad alignments encoded as strings
        final JavaRDD<AssemblyContigWithFineTunedAlignments> informativeContigs = informativeAndNotSo._1;
        final Tuple2<JavaRDD<AssemblyContigWithFineTunedAlignments>, JavaRDD<AssemblyContigWithFineTunedAlignments>> split =
                RDDUtils.split(informativeContigs, TempMultipleAlignmentsReclassifier::hasOnly2GoodAlignments, false);
        final JavaRDD<AssemblyContigWithFineTunedAlignments> twoGoodAlignments = split._1;

        // assembly contigs with more than 2 good alignments: without and with a complete picture
        final JavaRDD<AssemblyContigWithFineTunedAlignments> multipleAlignments = split._2;
        final Tuple2<JavaRDD<AssemblyContigWithFineTunedAlignments>, JavaRDD<AssemblyContigWithFineTunedAlignments>> split1 =
                RDDUtils.split(multipleAlignments, TempMultipleAlignmentsReclassifier::hasIncomePicture, false);
        final JavaRDD<AssemblyContigWithFineTunedAlignments> multipleAlignmentsIncompletePicture = split1._1;
        final JavaRDD<AssemblyContigWithFineTunedAlignments> multipleAlignmentsCompletePicture = split1._2;

        return new MultipleAlignmentReclassificationResults(garbage, multipleAlignmentsIncompletePicture,
                twoGoodAlignments, multipleAlignmentsCompletePicture);
    }

    // below are dummy predicates

    private static boolean isInformative(final AssemblyContigWithFineTunedAlignments contig) {
        return contig.isInformative();
    }

    private static boolean hasOnly2GoodAlignments(final AssemblyContigWithFineTunedAlignments contig) {
        return contig.updatedContig.alignmentIntervals.size() == 2;
    }

    private static boolean hasIncomePicture(final AssemblyContigWithFineTunedAlignments contig) {
        return contig.hasIncomePicture();
    }

    // =================================================================================================================

    /**
     * Process provided {@code originalConfiguration} of an assembly contig and split between good and bad alignments.
     * The returned pair has good alignments as its first (i.e. updated configuration),
     * and bad ones as its second (could be used for, e.g. annotate mappings of inserted sequence).
     *
     * <p>
     *     What is considered good and bad?
     *     For a particular mapping/alignment, it may offer low uniqueness in two sense:
     *     <ul>
     *         <li>
     *             low REFERENCE UNIQUENESS: meaning the sequence being mapped to match multiple locations on the reference;
     *         </li>
     *         <li>
     *             low READ UNIQUENESS: with only a very short part of the read being uniquely explained by this particular alignment;
     *         </li>
     *     </ul>
     *     Good alignments offer both high reference uniqueness and read uniqueness, as judged by the requested
     *     {@code mapQThresholdInclusive} and {@code uniqReadLenInclusive}
     *     (yes we are doing a hard-filtering but more advanced model is not our priority right now 2017-11-20).
     * </p>
     *
     * Note that "original" is meant to be possibly different from the returned configuration,
     * but DOES NOT mean the alignments of the contig as given by the aligner,
     * i.e. the configuration should be one of the best given by
     * {@link FilterLongReadAlignmentsSAMSpark#pickBestConfigurations(AlignedContig, Set, Double)}.
     */
    static Tuple2<List<AlignmentInterval>, List<AlignmentInterval>> removeNonUniqueMappings(final List<AlignmentInterval> originalConfiguration,
                                                                                            final int mapQThresholdInclusive,
                                                                                            final int uniqReadLenInclusive) {
        Utils.validateArg(originalConfiguration.size() > 2,
                "assumption that input configuration to be fine tuned has more than 2 alignments is violated.\n" +
                        originalConfiguration.stream().map(AlignmentInterval::toPackedString).collect(Collectors.toList()));

        // two pass, each focusing on removing the alignments of a contig that offers low uniqueness in one sense:

        // first pass is for removing alignments with low REFERENCE UNIQUENESS, using low mapping quality as the criterion
        final List<AlignmentInterval> selectedAlignments = new ArrayList<>(originalConfiguration.size()),
                                      lowUniquenessMappings = new ArrayList<>(originalConfiguration.size());

        for (final AlignmentInterval alignment : originalConfiguration) {
            if (alignment.mapQual >= mapQThresholdInclusive)
                selectedAlignments.add(alignment);
            else
                lowUniquenessMappings.add(alignment);
        }

        // second pass, the slower one, is to remove alignments offering low READ UNIQUENESS,
        // i.e. with only a very short part of the read being uniquely explained by this particular alignment;
        // the steps are:
        //      search bi-directionally until cannot find overlap any more,
        //      subtract the overlap from the distance covered on the contig by the alignment.
        //      This gives unique read region it explains.
        //      If this unique read region is "short": shorter than {@code uniqReadLenInclusive}), drop it.

        // each alignment has an entry of a tuple2, one for max overlap front, one for max overlap rear,
        // max overlap front is a tuple2 registering the index and overlap bases count
        final List<Tuple2<Tuple2<Integer, Integer>, Tuple2<Integer, Integer>>> maxOverlapMap =
                getMaxOverlapPairs(selectedAlignments);

        final List<Integer> idxToRemove = new ArrayList<>(selectedAlignments.size());
        for (int i = 0; i < selectedAlignments.size(); ++i) {
            final AlignmentInterval cur = selectedAlignments.get(i);
            final Tuple2<Tuple2<Integer, Integer>, Tuple2<Integer, Integer>> maxOverlapFrontAndRear = maxOverlapMap.get(i);
            final int maxOverlapFront = Math.max(0, maxOverlapFrontAndRear._1._2);
            final int maxOverlapRead = Math.max(0, maxOverlapFrontAndRear._2._2);

            final int uniqReadSpan = cur.endInAssembledContig - cur.startInAssembledContig + 1 - maxOverlapFront - maxOverlapRead;
            if (uniqReadSpan < uniqReadLenInclusive)
                idxToRemove.add(i);
        }

        if ( idxToRemove.isEmpty() )
            return new Tuple2<>(selectedAlignments, lowUniquenessMappings);

        // removing in reverse order so that iterators are not invalidated if we were to remove from start
        final ListIterator<Integer> rit = idxToRemove.listIterator(idxToRemove.size());
        while (rit.hasPrevious()) {
            selectedAlignments.remove( rit.previous().intValue() );
        }

        return new Tuple2<>(selectedAlignments, lowUniquenessMappings);
    }

    /**
     * Each alignment has an entry of a tuple2, one for max overlap front, one for max overlap rear,
     * max overlap front is a tuple2 registering the index and overlap bases count
     */
    private static List<Tuple2<Tuple2<Integer, Integer>, Tuple2<Integer, Integer>>> getMaxOverlapPairs(final List<AlignmentInterval> selectedAlignments) {

        final List<Tuple2<Tuple2<Integer, Integer>, Tuple2<Integer, Integer>>> maxOverlapMap =
                new ArrayList<>(Collections.nCopies(selectedAlignments.size(), new Tuple2<>(new Tuple2<>(-1, -1), new Tuple2<>(-1, -1))));

        for(int i = 0; i < selectedAlignments.size() - 1; ++i) {
            final AlignmentInterval cur = selectedAlignments.get(i);
            int maxOverlap = -1;
            int maxOverlapIdx = -1;
            for (int j = i + 1; j < selectedAlignments.size(); ++j) {
                final int overlap = AlignmentInterval.overlapOnContig(cur, selectedAlignments.get(j));
                if (overlap > 0) {
                    maxOverlap = Math.max(maxOverlap, overlap);
                    maxOverlapIdx = j;
                } else { // following ones, as guaranteed by the ordering of alignments in the contig, cannot overlap
                    break;
                }
            }
            if (maxOverlap > 0){
                // first set the max_overlap_rear of the current alignment
                Tuple2<Tuple2<Integer, Integer>, Tuple2<Integer, Integer>> oldValue = maxOverlapMap.get(i);
                maxOverlapMap.set(i, new Tuple2<>(oldValue._1, new Tuple2<>(maxOverlapIdx, maxOverlap))); // maxOverlapIdx cannot be -1 here
                // then conditionally set the max_overlap_front of the maxOverlapIdx-th alignment that maximally overlaps with the current alignment
                oldValue = maxOverlapMap.get(maxOverlapIdx);
                if (oldValue._1._2 < maxOverlap)
                    maxOverlapMap.set(maxOverlapIdx, new Tuple2<>(new Tuple2<>(i, maxOverlap), oldValue._2));
            }
        }
        return maxOverlapMap;
    }
}
