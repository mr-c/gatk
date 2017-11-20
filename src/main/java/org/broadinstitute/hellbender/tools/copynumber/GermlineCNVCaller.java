package org.broadinstitute.hellbender.tools.copynumber;

import org.apache.commons.collections4.ListUtils;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.SimpleCount;
import org.broadinstitute.hellbender.tools.copynumber.coverage.readcount.SimpleCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.RecordCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.io.Serializable;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Denoise and call copy-number variants in germline samples given their fragment counts and the corresponding output
 * of {@link DetermineGermlineContigPloidy}.
 * The former should be either HDF5 or TSV count files generated by {@link CollectFragmentCounts}.
 * Depending on available memory, it may be necessary to run over a subset of the intervals,
 * which can be specified by -L and must be present in all of the count files.
 *
 *  <p>
 *      If multiple samples are input, then the output is a model directory (which can be used for subsequently denoising
 *      individual samples, see below), as well as directories containing files that specify calls and
 *      sample-level model parameters for each sample.  The latter should be used as input to
 *      {@link //TODO PostprocessGermlineCNVCalls}.
 *  </p>
 *
 *  <p>
 *      If a single sample and a model directory are input, then only files that specify calls and sample-level
 *      model parameters for that sample are output.  Again, these should be used as input to
 *      {@link //TODO PostprocessGermlineCNVCalls}.
 *  </p>
 *
 * TODO Mehrtash can add documentation here
 *
 * <h3>Examples</h3>
 *
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" GermlineCNVCaller \
 *   -L intervals.interval_list \
 *   --input normal_1.counts.hdf5 \
 *   --input normal_2.counts.hdf5 \
 *   ... \
 *   --output output_dir \
 *   --outputPrefix normal_cohort
 * </pre>
 *
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" GermlineCNVCaller \
 *   -L intervals.interval_list \
 *   --input normal_1.counts.hdf5 \
 *   --model normal_cohort.ploidyModel.tsv \
 *   ... \
 *   --output output_dir \
 *   --outputPrefix normal_1
 * </pre>
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Denoise and call copy-number variants in germline samples given their fragment counts and the output of DetermineGermlineContigPloidy.",
        oneLineSummary = "Denoise and call copy-number variants in germline samples given their fragment counts and the output of DetermineGermlineContigPloidy.",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class GermlineCNVCaller extends CommandLineProgram {
    private enum Mode {
        COHORT, CASE
    }

    private static final String COHORT_DENOISING_CALLING_PYTHON_SCRIPT = "cohort_denoising_calling.py";
    private static final String CASE_SAMPLE_CALLING_PYTHON_SCRIPT = "case_sample_calling.py";

    private static final String OUTPUT_CALLS_SUFFIX = "-calls";

    public static final String READ_DEPTH_TABLE_FILE_LONG_NAME = "readDepth";
    public static final String READ_DEPTH_TABLE_FILE_SHORT_NAME = "depth";

    public static final String CONTIG_PLOIDY_TABLE_FILE_LONG_NAME = "contigPloidy";
    public static final String CONTIG_PLOIDY_TABLE_FILE_SHORT_NAME = "ploidy";

    @Argument(
            doc = "Input read-count files containing integer read counts in genomic intervals for all samples.  " +
                    "Intervals must be identical and in the same order for all samples.  " +
                    "If only a single sample is specified, a model directory must also be specified.  ",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            minElements = 1
    )
    private List<File> inputReadCountFiles = new ArrayList<>();

    @Argument(
            doc = "Input read-depth table file (output of DetermlineGermlineContigPloidy).",
            fullName = READ_DEPTH_TABLE_FILE_LONG_NAME,
            shortName = READ_DEPTH_TABLE_FILE_SHORT_NAME
    )
    private File inputReadDepthTableFile;

    @Argument(
            doc = "Input contig-ploidy table file (output of DetermlineGermlineContigPloidy).",
            fullName = CONTIG_PLOIDY_TABLE_FILE_LONG_NAME,
            shortName = CONTIG_PLOIDY_TABLE_FILE_SHORT_NAME
    )
    private File inputContigPloidyTableFile;

    @Argument(
            doc = "Input denoising-model directory.  If only a single sample is specified, this model will be used.  " +
                    "If multiple samples are specified, a new model will be built and this input will be ignored.",
            fullName = CopyNumberStandardArgument.MODEL_LONG_NAME,
            shortName = CopyNumberStandardArgument.MODEL_SHORT_NAME,
            optional = true
    )
    private String modelDir = null;

    @Argument(
            doc = "Input annotated-interval file containing annotations for GC content in genomic intervals (output of AnnotateIntervals).  " +
                    "All intervals specified via -L must be contained.  " +
                    "This input will be ignored in case-calling mode.",
            fullName = CopyNumberStandardArgument.ANNOTATED_INTERVALS_FILE_LONG_NAME,
            shortName = CopyNumberStandardArgument.ANNOTATED_INTERVALS_FILE_SHORT_NAME,
            optional = true
    )
    private File annotatedIntervalsFile = null;

    @Argument(
            doc = "Prefix for output filenames.",
            fullName =  CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME,
            shortName = CopyNumberStandardArgument.OUTPUT_PREFIX_SHORT_NAME
    )
    private String outputPrefix;

    @Argument(
            doc = "Output directory.",
            fullName =  StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private String outputDir;

    @Advanced
    @ArgumentCollection
    private GermlineDenoisingArgumentCollection germlineDenoisingArgumentCollection = new GermlineDenoisingArgumentCollection();

    @Advanced
    @ArgumentCollection
    private GermlineCallingArgumentCollection germlineCallingArgumentCollection = new GermlineCallingArgumentCollection();

    private Mode mode;

    @Override
    protected Object doWork() {
        setModeAndValidateArguments();

        //read in count files, validate they contain specified subset of intervals, and output count files for these intervals to temporary files
        final List<File> intervalSubsetReadCountFiles = writeIntervalSubsetReadCountFiles(inputReadCountFiles, )

        //call python inference code
        final boolean pythonReturnCode = executeGermlineCNVCallerPythonScript(intervalSubsetReadCountFiles, intervalsFile);

        if (!pythonReturnCode) {
            throw new UserException("Python return code was non-zero.");
        }

        logger.info("Germline denoising and CNV calling complete.");

        return "SUCCESS";
    }

    private void setModeAndValidateArguments() {
        if (modelDir == null) {
            if (inputReadCountFiles.size() > 1) {
                logger.info("Multiple samples provided, running in cohort mode...");
                mode = Mode.COHORT;
            } else {
                throw new UserException("Multiple samples must be provided if a denoising-model directory is not.");
            }
        } else {
            Utils.validateArg(!new File(modelDir).exists(), "Denoising-model directory does not exist.");
            if (inputReadCountFiles.size() > 1) {
                logger.warn("Multiple samples and a denoising-model directory were provided; the latter will be ignored...");
                mode = Mode.COHORT;
            } else {
                logger.info("A single sample and a denoising-model directory were provided, running in case mode...");
                mode = Mode.CASE;
            }
        }

        Utils.validateArg(inputReadCountFiles.size() == new HashSet<>(inputReadCountFiles).size(),
                "List of input read-count files cannot contain duplicates.");
        inputReadCountFiles.forEach(IOUtils::canReadFile);

        //TODO change this once output of DetermineGerlineContigPloidy is changed to directories
        IOUtils.canReadFile(inputReadDepthTableFile);
        IOUtils.canReadFile(inputContigPloidyTableFile);

        if (annotatedIntervalsFile != null) {
            IOUtils.canReadFile(annotatedIntervalsFile);
            if (mode == Mode.CASE) {
                logger.warn("Running in case mode, but an annotated-intervals file was provided; it will be ignored...");
            }
        }

        Utils.nonNull(outputPrefix);
        if (!new File(outputDir).exists()) {
            throw new UserException(String.format("Output directory %s does not exist.", outputDir));
        }

        //TODO validate argument collections
    }

    private boolean executeGermlineCNVCallerPythonScript(final List<File> intervalSubsetReadCountFiles,
                                                         final File intervalsFile) {
        final PythonScriptExecutor executor = new PythonScriptExecutor(true);
        final String outputDirArg = Utils.nonEmpty(outputDir).endsWith(File.separator) ? outputDir : outputDir + File.separator;    //add trailing slash if necessary
        final List<String> arguments = new ArrayList<>(Arrays.asList(
                "--modeling_interval_list=" + intervalsFile.getAbsolutePath(),
                "--sample_read_depth_metadata_table=" + inputReadDepthTableFile.getAbsolutePath(),
                "--sample_ploidy_metadata_table=" + inputContigPloidyTableFile.getAbsolutePath(),
                "--output_calls_path=" + outputDirArg + outputPrefix + OUTPUT_CALLS_SUFFIX));
        arguments.addAll(germlineDenoisingArgumentCollection.generatePythonArguments());
        arguments.addAll(germlineCallingArgumentCollection.generatePythonArguments());

        if (annotatedIntervalsFile != null) {
            arguments.add("--enable_explicit_gc_bias_modeling True");
        }

        if (mode == Mode.COHORT) {
            arguments.add("--output_model_path=" + outputDirArg + outputPrefix);
            return executor.executeScript(
                    new Resource(COHORT_DENOISING_CALLING_PYTHON_SCRIPT, GermlineCNVCaller.class),
                    null,
                    arguments);
        } else {
            return executor.executeScript(
                    new Resource(CASE_SAMPLE_CALLING_PYTHON_SCRIPT, GermlineCNVCaller.class),
                    null,
                    arguments);
        }
    }

    private static final class GermlineDenoisingArgumentCollection implements Serializable {
        private enum GCExpectationMode {
            MAP("map"),
            EXACT("exact"),
            HYBRID("hybrid");

            final String pythonArgumentString;

            GCExpectationMode(final String pythonArgumentString) {
                this.pythonArgumentString = pythonArgumentString;
            }
        }

        private static final long serialVersionUID = 1L;

        @Argument(
                doc = "Maximum number of bias factors.",
                fullName = "maxBiasFactors",
                minValue = 0,
                optional = true
        )
        private int maxBiasFactors = 5;

        @Argument(
                doc = "Typical mapping error rate.",
                fullName = "mappingErrorRate",
                minValue = 0.,
                optional = true
        )
        private double mappingErrorRate = 0.01;

        @Argument(
                doc = "Typical scale of interval-specific unexplained variance.",
                fullName = "intervalPsiScale",
                minValue = 0.,
                optional = true
        )
        private double intervalPsiScale = 0.001;

        @Argument(
                doc = "Typical scale of sample-specific unexplained variance.",
                fullName = "samplePsiScale",
                minValue = 0.,
                optional = true
        )
        private double samplePsiScale = 0.0001;

        @Argument(
                doc = "Precision of read depth pinning to its global value.",
                fullName = "depthCorrectionTau",
                minValue = 0.,
                optional = true
        )
        private double depthCorrectionTau = 10000.0;

        @Argument(
                doc = "Standard deviation of log mean bias.",
                fullName = "logMeanBiasStandardDeviation",
                minValue = 0.,
                optional = true
        )
        private double logMeanBiasStandardDeviation = 0.1;

        @Argument(
                doc = "Initial value of ARD prior precision relative to the typical interval-specific unexplained variance scale.",
                fullName = "initARDRelUnexplainedVariance",
                minValue = 0.,
                optional = true
        )
        private double initARDRelUnexplainedVariance = 0.1;

        @Argument(
                doc = "Number of knobs on the GC curves.",
                fullName = "numGCBins",
                minValue = 1,
                optional = true
        )
        private int numGCBins = 20;

        @Argument(
                doc = "Prior standard deviation of the GC curve from flat.",
                fullName = "gcCurveStandardDeviation",
                minValue = 0.,
                optional = true
        )
        private double gcCurveStandardDeviation = 1.;

        @Argument(
                doc = "The strategy for calculating copy number posterior expectations in the denoising model.",
                fullName = "gcExpectationMode",
                optional = true
        )
        private GCExpectationMode gcExpectationMode = GCExpectationMode.HYBRID;

        @Argument(
                doc = "Enable discovery of bias factors.",
                fullName = "enableBiasFactors",
                optional = true
        )
        private boolean enableBiasFactors = true;

        @Argument(
                doc = "Disable bias factor discovery in intervals in CNV-active regions.",
                fullName = "disableBiasFactorsInFlatClass",
                optional = true
        )
        private boolean disableBiasFactorsInFlatClass = false;

        private List<String> generatePythonArguments() {
            return Arrays.asList(
                    String.format("--max_bias_factors=%d", maxBiasFactors),
                    String.format("--mapping_error_rate=%f", mappingErrorRate),
                    String.format("--psi_i_scale=%f", intervalPsiScale),
                    String.format("--psi_s_scale=%f", samplePsiScale),
                    String.format("--depth_correction_tau=%f", depthCorrectionTau),
                    String.format("--log_mean_bias_std=%f", logMeanBiasStandardDeviation),
                    String.format("--init_ard_rel_unexplained_variance=%f", initARDRelUnexplainedVariance),
                    String.format("--num_gc_bins=%d", numGCBins),
                    String.format("--gc_curve_sd=%f", gcCurveStandardDeviation),
                    String.format("--q_c_expectation_mode=%s", gcExpectationMode.pythonArgumentString),
                    String.format("--enable_bias_factors=%s", enableBiasFactors ? "True" : "False"),
                    String.format("--disable_bias_factors_in_flat_class=%s", disableBiasFactorsInFlatClass ? "True" : "False"));
        }
    }

    private static final class GermlineCallingArgumentCollection implements Serializable {
        private static final long serialVersionUID = 1L;

        @Argument(
                doc = "Prior probability of alt copy number with respect to contig baseline state in the reference copy number.",
                fullName = "pAlt",
                minValue = 0.,
                optional = true
        )
        private double pAlt = 1E-6;

        @Argument(
                doc = "Prior probability of using flat copy number distribution as prior.",
                fullName = "pFlat",
                minValue = 0.,
                optional = true
        )
        private double pFlat = 1E-3;

        @Argument(
                doc = "Coherence length of CNV events (in the units of bp).",
                fullName = "cnvCoherenceLength",
                minValue = 0.,
                optional = true
        )
        private double cnvCoherenceLength = 10000.0;

        @Argument(
                doc = "Coherence length of copy number classes (in the units of bp).",
                fullName = "classCoherenceLength",
                minValue = 0.,
                optional = true
        )
        private double classCoherenceLength = 10000.0;

        @Argument(
                doc = "Highest considered copy number.",
                fullName = "maxCopyNumber",
                minValue = 0,
                optional = true
        )
        private int maxCopyNumber = 5;

        @Argument(
                doc = "Initialize with flat copy number prior everywhere.",
                fullName = "doInitializeToFlatClass",
                optional = true
        )
        private boolean doInitializeToFlatClass = true;

        private List<String> generatePythonArguments() {
            return Arrays.asList(
                    String.format("--p_alt=%f", pAlt),
                    String.format("--p_flat=%f", pFlat),
                    String.format("--cnv_coherence_length=%f", cnvCoherenceLength),
                    String.format("--class_coherence_length=%f", classCoherenceLength),
                    String.format("--max_copy_number=%d", maxCopyNumber),
                    String.format("--initialize_to_flat_class=%s", doInitializeToFlatClass ? "True" : "False"));
        }
    }
}
