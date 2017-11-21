package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.Serializer;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.ChimericAlignment;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.StrandSwitch;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

@DefaultSerializer(Serializer.class)
final class AssemblyContigWithFineTunedAlignments {
    private static final AlignedContig.Serializer contigSerializer = new AlignedContig.Serializer();
    public static final List<String> emptyInsertionMappings = Collections.emptyList();

    final AlignedContig contig;
    final List<String> insertionMappings;

    AssemblyContigWithFineTunedAlignments(final AlignedContig contig) {
        this(contig, emptyInsertionMappings);
    }

    AssemblyContigWithFineTunedAlignments(final AlignedContig contig,
                                          final List<String> insertionMappings) {
        this.contig = contig;
        this.insertionMappings = Utils.nonNull(insertionMappings);
    }

    AssemblyContigWithFineTunedAlignments(final Kryo kryo, final Input input) {
        contig = contigSerializer.read(kryo, input, AlignedContig.class);
        final int insertionMappingsSize = input.readInt();
        insertionMappings = new ArrayList<>(insertionMappingsSize);
        if (insertionMappingsSize != 0) {
            for (int i = 0; i < insertionMappingsSize; ++i) {
                insertionMappings.add(input.readString());
            }
        }
    }

    // TODO: 11/21/17 push these predicates to AlignedContig class in a later PR
    // after fine tuning, a contig may have no good alignment left, or only 1
    final boolean isInformative() {
        return contig.alignmentIntervals.size() > 1;
    }

    final boolean hasOnly2GoodAlignments () {
        return contig.alignmentIntervals.size() == 2;
    }

    final boolean hasIncomePicture() {
        if ( hasOnly2GoodAlignments() )
            return hasIncompletePictureFromTwoAlignments();
        else
            return hasIncompletePictureFromMultipleAlignments();
    }

    //==================================================================================================================

    boolean hasIncompletePictureFromTwoAlignments() {
        return hasIncompletePictureDueToRefSpanContainment()
                ||
                (firstAndLastAlignmentMappedToSameChr() && hasIncompletePictureDueToOverlappingRefOrderSwitch());
    }

    boolean hasIncompletePictureDueToRefSpanContainment() {

        final AlignmentInterval one = contig.alignmentIntervals.get(0);
        final AlignmentInterval two = contig.alignmentIntervals.get(1);
        return one.referenceSpan.contains(two.referenceSpan)
                ||
                two.referenceSpan.contains(one.referenceSpan);
    }

    boolean hasIncompletePictureDueToOverlappingRefOrderSwitch() {

        final AlignmentInterval one = contig.alignmentIntervals.get(0);
        final AlignmentInterval two = contig.alignmentIntervals.get(1);
        final SimpleInterval referenceSpanOne = one.referenceSpan;
        final SimpleInterval referenceSpanTwo = two.referenceSpan;

        if (referenceSpanOne.contains(referenceSpanTwo) || referenceSpanTwo.contains(referenceSpanOne))
            return true;

        // TODO: 10/29/17 this obsoletes the inverted duplication call code we have now, but those could be used to figure out how to annotate which known ref regions are invert duplicated
        if (one.forwardStrand != two.forwardStrand) {
            return referenceSpanOne.overlaps(referenceSpanTwo);
        } else {
            if (one.forwardStrand) {
                return referenceSpanOne.getStart() > referenceSpanTwo.getStart() &&
                        referenceSpanOne.getStart() <= referenceSpanTwo.getEnd();
            } else {
                return referenceSpanTwo.getStart() > referenceSpanOne.getStart() &&
                        referenceSpanTwo.getStart() <= referenceSpanOne.getEnd();
            }
        }
    }

    boolean firstAndLastAlignmentMappedToSameChr() {
        Utils.validateArg(!hasIncompletePictureDueToRefSpanContainment(),
                "assumption that input contig has alignments whose ref span is not completely covered by the other's is violated \n" +
                        contig.toString());

        final String firstMappedChr = this.contig.alignmentIntervals.get(0).referenceSpan.getContig();
        final String lastMappedChr  = this.contig.alignmentIntervals.get(this.contig.alignmentIntervals.size() - 1).referenceSpan.getContig();

        return firstMappedChr.equals(lastMappedChr);
    }

    //==================================================================================================================

    private boolean hasIncompletePictureFromMultipleAlignments() {
        final int goodAlignmentsCount = contig.alignmentIntervals.size();
        Utils.validateArg(goodAlignmentsCount > 2,
                "assumption that input contig has more than 2 alignments is violated.\n" +
                        contig.toString());

        final AlignmentInterval head = contig.alignmentIntervals.get(0),
                                tail = contig.alignmentIntervals.get(goodAlignmentsCount-1);

        // tail not resuming the "direction of flow" started by head
        if ( !head.referenceSpan.getContig().equals(tail.referenceSpan.getContig()) || head.forwardStrand != tail.forwardStrand)
            return true;

        // reference span anomaly
        final SimpleInterval referenceSpanHead = head.referenceSpan,
                             referenceSpanTail = tail.referenceSpan;
        if (referenceSpanHead.contains(referenceSpanTail) || referenceSpanTail.contains(referenceSpanHead))
            return true;

        final boolean refSpanAnomaly =
                contig.alignmentIntervals
                        .subList(1, goodAlignmentsCount - 1)
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
        contigSerializer.write(kryo, output, contig);
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
