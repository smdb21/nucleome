package edu.scripps.yates.nucleome.turboID;

import static org.junit.Assert.fail;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.fitting.GaussianCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinRetriever;
import edu.scripps.yates.excel.ExcelReader;
import edu.scripps.yates.nucleome.turboID.annotations.AnnotationsUtil;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.maths.Histogram;
import edu.scripps.yates.utilities.maths.Maths;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import gnu.trove.set.hash.THashSet;

public class TurboIDDataAnalysisOptimalParams {
	private final static Logger log = Logger.getLogger(TurboIDDataAnalysisOptimalParams.class);
	private static final String ACCESSION = "accession";
	private static final String TRANSMEMBRANE = "Transmembrane";
	private static final String DESCRIPTION = "description";
	private static final String SPC = "spec count";
	private static final String PSM = "psm";
	private static final int GO_TERM_DISTANCE = 3;
	// input data with all the experiments in one excel file
	private final File excelFile_TMT11_TurboID_NuCy = new File(
			"Z:\\share\\Salva\\data\\4D_Nucleome\\TurboID\\input\\optimized_params\\counts_NE4_NuCy_reload_sch11str_noP.xlsx");
	// input data with all the experiments in one excel file with the difference of
	// having a filter for at least 2 unique peptides per protein
	private final File excelFile_TMT11_TurboID_NuCy_2unique = new File(
			"Z:\\share\\Salva\\data\\4D_Nucleome\\TurboID\\input\\optimized_params\\counts_NE4_NuCy_reload_sch11str_noP_2unique.xlsx");
	// TEMPORALY DISABLED
	private final String[] sheetsTMT11Nu = { "Nu_Ta_reload_sch11str", "Nu_Tb_reload_sch11str",
			"Nu_Tb_rerun_reload_sch11str" };
//	, "20191016_1114LT-Nu_Ta11_rerun" };
//	private final String[] sheetsTMT11Nu = { "20191016_1114LT-Nu_Ta11_rerun" };

	private final String[] sheetsTMT11Cy = { "Cy_Ta_reload_sch11str", "Cy_Tb_reload_sch11str" };

	private final File outputFolder = excelFile_TMT11_TurboID_NuCy.getParentFile();
	private final File uniprotReleasesFolder = new File("Z:\\share\\Salva\\data\\uniprotKB");
	private final AnnotationsUtil annotationsUtil = new AnnotationsUtil(uniprotReleasesFolder);
	private final boolean[] yesno = { true, false };
	// private boolean normalizedByMixChannel;
	// private boolean normalizeByV5Tag;
	private final THashMap<String, Replicate> replicateBySheets = new THashMap<String, Replicate>();
	private boolean normalizeChannelsManually;
	private double timesSigmaForSpecificity;
	private final DecimalFormat formatter = new DecimalFormat("#.#");

	// *************************
	//
	// GENERATE OUTPUT FOR SAINTExpress
	//
	private final boolean generateSAINTExpressOutput = false;
	private final boolean useJustNucleusFractionForSAINTExpress = true;
	private final boolean runAllBaitsTogether = true;
	// *************************

	// *************************
	// FRACTION TO ANALYZE
	// *************************
	private final boolean analyzeNucleusFraction = true;
	private final boolean analyzeCytoplasmFraction = !useJustNucleusFractionForSAINTExpress;// false;// true;

	// *************************
	// MINIMUM THRESHOLD FOR RAW INTENSITIES
	private static final Double minimumIntensityThreshold = null;// 5000.0;
	// *************************
	//////////////////////////////

	public static void main(String[] args) {
		try {
			new TurboIDDataAnalysisOptimalParams().run();
		} catch (final Exception e) {
			e.printStackTrace();
			System.exit(-1);
		}
		System.exit(0);
	}

	public void run() {
		try {
			if (!useJustNucleusFractionForSAINTExpress && generateSAINTExpressOutput && !analyzeCytoplasmFraction) {
				throw new IllegalArgumentException(
						"You have to enable analyzeCytoplasmFraction! in order to generate SAINTExpress output");
			}
			if (generateSAINTExpressOutput && !analyzeNucleusFraction) {
				throw new IllegalArgumentException(
						"You have to enable analyzeNucleusFraction! in order to generate SAINTExpress output");
			}
			setupChannels();
			final boolean applyLog2 = true;

			// for (final boolean useNormalizedIntensities : yesno) {
			final boolean onlyTM = false;
			final boolean onlyNonTM = false;
			final boolean onlyComplete = true; // only use complete proteins, which means that are present in 4 out of 6
												// replicates.
			final boolean useNormalizedIntensities = true;
			final boolean ignoreMissingValues = true;
			// normalizations:
			final boolean normalizedByMixChannel = false; // discarded, always
															// FALSE
			final boolean normalizeByV5Tag = false; // discarded, always FALSE
			final boolean normalizeToAverage = false;
			final boolean normalizeToTurboIDChannelAverage = false;
			normalizeChannelsManually = true; // normalize by the sum of the
												// channel
			timesSigmaForSpecificity = 2;

			String fileoutputPrefix = null;
			int minPerGroup = 0; // minimum number of replicates that is required to be a protein detected (or
									// also, per bait, depending on where this param is used)

			final List<TurboIDExperiment> experiments = new ArrayList<TurboIDExperiment>();
			if (analyzeNucleusFraction) {
				final TurboIDExperiment turboIDExperiment = readTurboIDExperiment(excelFile_TMT11_TurboID_NuCy_2unique,
						sheetsTMT11Nu, TurboIDFraction.NU);
				fileoutputPrefix = "Nu";
				minPerGroup = 4;
				experiments.add(turboIDExperiment);
			}
			if (analyzeCytoplasmFraction) {
				final TurboIDExperiment turboIDExperiment = readTurboIDExperiment(excelFile_TMT11_TurboID_NuCy_2unique,
						sheetsTMT11Cy, TurboIDFraction.CY);
				fileoutputPrefix = "Cy";
				minPerGroup = 4;
				experiments.add(turboIDExperiment);
			}
			final File outputFile = new File(outputFolder.getAbsolutePath() + File.separator
					+ getOutputFileName(fileoutputPrefix, onlyTM, onlyNonTM, onlyComplete, ignoreMissingValues,
							useNormalizedIntensities, normalizeToAverage, normalizeToTurboIDChannelAverage, applyLog2));
			for (final TurboIDExperiment turboIDExperiment : experiments) {

				final List<String> accs = new ArrayList<String>();
				turboIDExperiment.getProteinsSorted(TurboIDExperiment.getComparatorByTotalSPC(true))
						.forEach(protein -> accs.add(protein.getAcc()));
				loadUniprotAnnotations(accs);
				System.out.println("");
				// NORMALIZE by mix
				if (normalizedByMixChannel) {
					turboIDExperiment.normalizeByMixChannel();
				}
				if (normalizeToAverage) {
					turboIDExperiment.normalizeNormalizedByAvgOfBaitChannels(ignoreMissingValues, applyLog2);
				}
				if (normalizeToTurboIDChannelAverage) {
					final Collection<Replicate> replicatesToUseOnNormalization = new THashSet<Replicate>();
					replicatesToUseOnNormalization.add(Replicate.A);
					turboIDExperiment.normalizeToTurboIDChannelAverage(ignoreMissingValues, applyLog2,
							replicatesToUseOnNormalization);
				}
				// NORMALIZE by V5 tag
				if (normalizeByV5Tag) {
					turboIDExperiment.normalizeByV5Tag(applyLog2);
				}

				// }
				// }
				// }
				// }
				if (!generateSAINTExpressOutput) {
					final double sigmaFactor = 0.25;
					for (final Replicate replicate : Replicate.values(turboIDExperiment.getFraction())) {
						final File fileGaussianRep = new File(outputFile.getParentFile().getAbsoluteFile()
								+ File.separator + FilenameUtils.getBaseName(outputFile.getAbsolutePath())
								+ "_GaussFit_" + replicate.name() + ".txt");
						calculateAndFitNoise(turboIDExperiment, replicate, fileGaussianRep, sigmaFactor);
					}
					final File fileGaussian = new File(outputFile.getParentFile().getAbsoluteFile() + File.separator
							+ FilenameUtils.getBaseName(outputFile.getAbsolutePath()) + "_GaussFit.txt");
					calculateAndFitNoise(turboIDExperiment, fileGaussian, sigmaFactor);

					final File outputFile5 = new File(outputFile.getParentFile().getAbsoluteFile() + File.separator
							+ FilenameUtils.getBaseName(outputFile.getAbsolutePath()) + "_specific.txt");
					exportSpecificity(turboIDExperiment, outputFile5);

					turboIDExperiment.exportToFileForClustering(outputFile, onlyTM, onlyNonTM, ignoreMissingValues,
							useNormalizedIntensities);

					// calculate ANOVA per protein comparing the baits

					final double anovaPValueThreshold = 0.05;

					// now the same but with distribSPC
//					final File outputFile8 = new File(outputFile.getParentFile().getAbsoluteFile() + File.separator
//							+ fileoutputPrefix + "_anova_log_distribSPC.txt");
//					calculateAnovaPerProteinComparingBaits(turboIDExperiment, outputFile8, true, true, false,
//							anovaPValueThreshold, minPerGroup);
					// now the same but with log distribInt
//					final File outputFile9 = new File(outputFile.getParentFile().getAbsoluteFile() + File.separator
//							+ fileoutputPrefix + "_anova_log_distribInt.txt");
//					calculateAnovaPerProteinComparingBaits(turboIDExperiment, outputFile9, true, false, true,
//							anovaPValueThreshold, minPerGroup);
					final File outputFile10 = new File(outputFile.getParentFile().getAbsoluteFile() + File.separator
							+ fileoutputPrefix + "_anova_distribInt.txt");
					calculateAnovaPerProteinComparingBaits(turboIDExperiment, outputFile10, false, false, true,
							anovaPValueThreshold, minPerGroup);
				}
			}
			if (generateSAINTExpressOutput) {
				final File outputFolderForSAINTExpress = new File(
						outputFile.getParentFile().getAbsolutePath() + File.separator + "SAINTExpress");
				final File outputFolderForUMAP = new File(
						outputFile.getParentFile().getAbsolutePath() + File.separator + "UMAP");
				TurboIDExperiment control = null;
				TurboIDExperiment test = null;
				for (final TurboIDExperiment turboIDExperiment2 : experiments) {
					if (turboIDExperiment2.getFraction() == TurboIDFraction.CY) {
						control = turboIDExperiment2;
					} else if (turboIDExperiment2.getFraction() == TurboIDFraction.NU) {
						test = turboIDExperiment2;
					}
				}
				final SAINTExpress saintExpress = new SAINTExpress(outputFolderForSAINTExpress, control, test);
				saintExpress.setRunAllBaitsTogether(runAllBaitsTogether);
				saintExpress.run(onlyTM, onlyNonTM, onlyComplete, minPerGroup);
				try {
					final UMAPClustering umap = new UMAPClustering(outputFolderForUMAP, control, test);
					umap.run(onlyTM, onlyNonTM);
				} catch (final IllegalArgumentException e) {
					log.warn(e.getMessage());
				}
			}
		} catch (final Exception e) {
			e.printStackTrace();
			fail();
		}

	}

	private void calculateAnovaPerProteinComparingBaits(TurboIDExperiment turboIDExperiment, File outputFile,
			boolean applyLog, boolean distribSPC, boolean distribIntensity, double anovaPValueThreshold,
			int minPerGroup) throws IOException {
		final File benchmarksFile = new File(
				outputFile.getParentFile().getAbsolutePath() + File.separator + "benchmarks_20191104.txt");
		final Map<String, String> benchmarks = Benchmarks.getBenchmarks(benchmarksFile);
		String prefix = "";
		if (distribSPC) {
			prefix = "distSPC";
		} else if (distribIntensity) {
			prefix = "distNormInt";
		} else {
			prefix = "normInt";
		}
		if (applyLog) {
			prefix = "log2_" + prefix;
		}
		log.info(prefix);
		final FileWriter fw = new FileWriter(outputFile);
		final TurboIDFraction fraction = turboIDExperiment.getFraction();
		final List<ProteinFromTurboID> proteins = turboIDExperiment.getProteinsSorted(
				TurboIDExperiment.getComparatorByAnova(true, applyLog, distribSPC, distribIntensity, fraction));
		fw.write("accession" + "\t" + "gene" + "\t" + "description" + "\t" + "transmembrane" + "\t" + "nucleus" + "\t"
				+ "DNA binding" + "\t" + "transcription factor" + "\t" + "RNA binding" + "\t" + "Heterochomatin" + "\t"
				+ "total SPCs" + "\t");
		// spc per replicate
		for (final Replicate replicate : Replicate.values(fraction)) {
			fw.write("SPC_" + replicate.name() + "\t");
		}
		fw.write("is complete" + "\t");
		fw.write("complete_4_6" + "\t");
		fw.write("replicates" + "\t");
		fw.write("benchmark" + "\t");
		fw.write(prefix + " ANOVA p-value" + "\t");
		fw.write(prefix + " ANOVA p-value (with TbID)" + "\t");
		// comparisons with t-test
		for (int i = 0; i < TurboIDExperimentType.getBaitsAndTBID().size(); i++) {
			final TurboIDExperimentType bait1 = TurboIDExperimentType.getBaitsAndTBID().get(i);
			for (int j = i + 1; j < TurboIDExperimentType.getBaitsAndTBID().size(); j++) {
				final TurboIDExperimentType bait2 = TurboIDExperimentType.getBaitsAndTBID().get(j);
				if (bait1 != bait2) {
					fw.write(bait1.name() + "-" + bait2.name() + "\t");
				}
			}
		}
		// distributed intensities averaged per bait
		for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaitsAndTBID()) {
			fw.write(prefix + " " + bait.name() + "\t");
		}
		// distributed intensities per channel
		for (final TurboIDExperimentType bait : TurboIDExperimentType.values()) {
			for (final Replicate replicate : Replicate.values(fraction)) {
				for (final TurboID_Channel_Norm channel : TurboID_Channel_Norm.values()) {
					if (channel.getExpType() == bait && channel.getReplicate() == replicate) {
						fw.write("distr_" + bait.name() + "_" + channel.name() + "\t");
					}
				}
			}
		}
		// normalized intensities per channel
		for (final TurboIDExperimentType bait : TurboIDExperimentType.values()) {
			for (final Replicate replicate : Replicate.values(fraction)) {
				for (final TurboID_Channel_Norm channel : TurboID_Channel_Norm.values()) {
					if (channel.getExpType() == bait && channel.getReplicate() == replicate) {
						fw.write("norm_" + bait.name() + "_" + channel.name() + "\t");
					}
				}
			}
		}
		// raw intensities per channel
		for (final TurboIDExperimentType bait : TurboIDExperimentType.values()) {
			for (final Replicate replicate : Replicate.values(fraction)) {
				for (final TurboID_Channel_Ori channel : TurboID_Channel_Ori.values()) {
					if (channel.getExpType() == bait && channel.getReplicate() == replicate) {
						fw.write("raw_" + bait.name() + "_" + channel.name() + "\t");
					}
				}
			}
		}
		// pseudo SPC per channel
		for (final TurboIDExperimentType bait : TurboIDExperimentType.values()) {
			for (final Replicate replicate : Replicate.values(fraction)) {
				for (final TurboID_Channel_Norm channel : TurboID_Channel_Norm.values()) {
					if (channel.getExpType() == bait && channel.getReplicate() == replicate) {
						fw.write("pseudoSPC_" + bait.name() + "_" + channel.name() + "\t");
					}
				}
			}
		}
		fw.write("\n");
		for (final ProteinFromTurboID protein : proteins) {
			if (protein.getAcc().equals("Q7TPN9")) {
				log.info(protein);
			}
			if (protein.getAcc().equals("Q8CHH9")) {
				System.out.println("äsdf");
			}
			fw.write(protein.getAcc() + "\t" + parseForExcelConflictingGenes(protein.getGene()) + "\t"
					+ annotationsUtil.getDescription(protein) + "\t" + annotationsUtil.getTransmembraneRegion(protein)
					+ "\t" + annotationsUtil.isNucleus(protein) + "\t" + annotationsUtil.isDNABinding(protein) + "\t"
					+ annotationsUtil.isTranscriptionFactor(protein) + "\t" + annotationsUtil.isRNABinding(protein)
					+ "\t" + annotationsUtil.isHeterochromatin(protein) + "\t" + protein.getSumSPCAcrossReplicates()
					+ "\t");
			// spc per replicate
			for (final Replicate replicate : Replicate.values(fraction)) {
				fw.write(protein.getSpc(replicate) + "\t");
			}
			// column saying whether contains enough data to do the pvalue or
			// not

			// is complete: 4/6 signals in all 4 baits
			final boolean isvalid = protein.isAnovaValid(distribSPC, distribIntensity, fraction, minPerGroup);
			fw.write(isvalid + "\t");
			// is complete_4_6: present in 4/6 replicates
			final boolean isvalid2 = protein.isComplete(minPerGroup, fraction);
			fw.write(isvalid2 + "\t");
			// num replicates
			final int numReps = protein.getNumReplicates(fraction);
			fw.write(numReps + "\t");
			// is benchmark
			final boolean isBenchmark = benchmarks.containsKey(protein.getGene());
			fw.write(isBenchmark + "\t");
			// ANOVA over the 4 baits
			final double anovaPValue = protein.getAnovaPValueOverBaits(applyLog, distribSPC, distribIntensity,
					fraction);
			fw.write(anovaPValue + "\t");
			// ANOVA over the 4 baits + turboID
			final double anovaPValue2 = protein.getAnovaPValueOverBaitsAndTurboIDOnly(applyLog, distribSPC,
					distribIntensity, fraction);
			fw.write(anovaPValue2 + "\t");
			// pairwise comparisons with t-test + TurboID
			for (int i = 0; i < TurboIDExperimentType.getBaitsAndTBID().size(); i++) {
				final TurboIDExperimentType bait1 = TurboIDExperimentType.getBaitsAndTBID().get(i);
				for (int j = i + 1; j < TurboIDExperimentType.getBaitsAndTBID().size(); j++) {
					final TurboIDExperimentType bait2 = TurboIDExperimentType.getBaitsAndTBID().get(j);
					if (bait1 != bait2) {
						final double ttestPValue = protein.getTTestPValueOverBaits(bait1, bait2, applyLog, distribSPC,
								distribIntensity, fraction);
						fw.write(ttestPValue + "\t");
					}
				}
			}
			// distributed intensities
			final TObjectDoubleHashMap<TurboID_Channel_Norm> distributedIntensitiesWithNormalizedIntensities = protein
					.getDistributedIntensitiesWithNormalizedIntensities();
			for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaitsAndTBID()) {
				double value = Double.NaN;
				if (distribSPC) {
					final TDoubleArrayList toAverage = new TDoubleArrayList();
					for (final Replicate replicate : Replicate.values(fraction)) {
						final TObjectDoubleHashMap<TurboID_Channel_Norm> pseudoSpecCountsWithNormalizedIntensities = protein
								.getPseudoSpecCountsWithNormalizedIntensities(replicate);
						for (final TurboID_Channel_Norm channel : pseudoSpecCountsWithNormalizedIntensities.keySet()) {
							if (channel.getExpType() == bait) {
								final double val = pseudoSpecCountsWithNormalizedIntensities.get(channel);
								if (!Double.isNaN(val)) {
									toAverage.add(val);
								}
							}
						}
					}
					value = Maths.mean(toAverage);
				} else if (distribIntensity) {
					final TDoubleArrayList toAverage = new TDoubleArrayList();

					for (final TurboID_Channel_Norm channel : distributedIntensitiesWithNormalizedIntensities
							.keySet()) {
						if (channel.getExpType() == bait) {
							final double val = distributedIntensitiesWithNormalizedIntensities.get(channel);
							if (!Double.isNaN(val)) {
								toAverage.add(val);
							}
						}
					}

					if (toAverage.isEmpty()) {
						value = 0.0;// Double.NaN;
					} else {
						value = Maths.mean(toAverage);
					}
				} else {
					value = protein.getAverageFromNormalized(bait, true);
				}
				if (applyLog) {
					fw.write(Maths.log(value, 2) + "\t");
				} else {
					fw.write(value + "\t");
				}
			}
			// print the distributed intensities
			for (final TurboIDExperimentType bait : TurboIDExperimentType.values()) {

				for (final Replicate replicate : Replicate.values(fraction)) {
					for (final TurboID_Channel_Norm channel : TurboID_Channel_Norm.values()) {
						if (channel.getExpType() == bait && channel.getReplicate() == replicate) {
							final double value = protein.getDistributedIntensitiesWithNormalizedIntensities()
									.get(channel);
							if (applyLog && !Double.isNaN(value)) {
								final double log2 = Maths.log(value, 2);
								fw.write(log2 + "\t");
							} else {
								fw.write(value + "\t");
							}

						}
					}
				}
			}
			// print the normalized intensities
			for (final TurboIDExperimentType bait : TurboIDExperimentType.values()) {

				for (final Replicate replicate : Replicate.values(fraction)) {
					for (final TurboID_Channel_Norm channel : TurboID_Channel_Norm.values()) {
						if (channel.getExpType() == bait && channel.getReplicate() == replicate) {
							final double value = protein.getNormalizedIntensities(bait).get(channel);
							if (applyLog && !Double.isNaN(value)) {
								final double log2 = Maths.log(value, 2);
								fw.write(log2 + "\t");
							} else {
								fw.write(value + "\t");
							}

						}
					}
				}
			}
			// print the raw intensities
			for (final TurboIDExperimentType bait : TurboIDExperimentType.values()) {
				for (final Replicate replicate : Replicate.values(fraction)) {
					for (final TurboID_Channel_Ori channel : TurboID_Channel_Ori.values()) {
						if (channel.getExpType() == bait && channel.getReplicate() == replicate) {
							final double value = protein.getOriginalIntensities(bait).get(channel);
							if (applyLog && !Double.isNaN(value)) {
								final double log2 = Maths.log(value, 2);
								fw.write(log2 + "\t");
							} else {
								fw.write(value + "\t");
							}

						}
					}
				}
			}
			// print the pseudoSPC
			for (final TurboIDExperimentType bait : TurboIDExperimentType.values()) {
				for (final Replicate replicate : Replicate.values(fraction)) {
					for (final TurboID_Channel_Norm channel : TurboID_Channel_Norm.values()) {
						if (channel.getExpType() == bait && channel.getReplicate() == replicate) {
							final double value = protein.getPseudoSpecCountsWithNormalizedIntensities(replicate)
									.get(channel);
							if (applyLog && !Double.isNaN(value)) {
								final double log2 = Maths.log(value, 2);
								fw.write(log2 + "\t");
							} else {
								fw.write(value + "\t");
							}

						}
					}
				}
			}
			fw.write("\n");
			if (anovaPValue >= anovaPValueThreshold) {
				// break;
			}
		}
		fw.close();
		log.info("File written at: " + outputFile.getAbsolutePath());

	}

	private void calculateAndFitNoise(TurboIDExperiment experiment, Replicate replicate, File outputFileGaussian,
			double sigmaFactorToFit) throws IOException {

		final List<WeightedObservedPoint> points = getHistogramOfLog2PseudoSPC(experiment, replicate);
		if (points.isEmpty()) {
			return;
		}
		experiment.setExperimentalHistogram(points, replicate);
		final double[] fit = GaussianCurveFitter.create().fit(points);
		final double normalization = fit[0];
		final double mean = fit[1];
		final double sigma = fit[2];
		System.out.println(normalization + "\t" + mean + "\t" + sigma);

		log.info("Now, I mirror the right part of the distribution from mean to max and then, I do fitting again");
		// final List<WeightedObservedPoint> newPoints =
		// getNewPointsToFitMirroringAfterInflexionPoint(points);
		log.info("Now, I take only up to " + sigmaFactorToFit + "x sigma to fit again the data");
		final List<WeightedObservedPoint> newPoints = getNewPointsToFitUpTo1Sigma(points, mean,
				sigmaFactorToFit * sigma);

		final double[] fit2 = GaussianCurveFitter.create().fit(newPoints);
		final double normalization2 = fit2[0];
		final double mean2 = fit2[1];
		final double sigma2 = fit2[2];
		System.out.println(normalization2 + "\t" + mean2 + "\t" + sigma2);

		final Gaussian gaussian = new Gaussian(normalization2, mean2, sigma2);
		final FileWriter fw = new FileWriter(outputFileGaussian);
		fw.write("Gaussian fitting:\n");
		fw.write("Normalization=" + normalization2 + "\n");
		fw.write("Mean=" + mean2 + "\n");
		fw.write("Sigma=" + sigma2 + "\n\n");
		fw.write("log2 pseudoSPC" + "\t" + "observed" + "\t" + "used for fitting" + "\t" + "fitted" + "\n");

		for (int index = 0; index < points.size(); index++) {
			final WeightedObservedPoint point = points.get(index);

			final double x = point.getX();
			final double observed = point.getY();
			String observedMirroredString = "";
			if (newPoints.size() > index) {
				final WeightedObservedPoint point2 = newPoints.get(index);
				observedMirroredString = String.valueOf(point2.getY());
			}
			final double fitted = gaussian.value(x);
			fw.write(x + "\t" + observed + "\t" + observedMirroredString + "\t" + fitted + "\n");
		}
		fw.close();
		log.info("File written at: " + outputFileGaussian.getAbsolutePath());

		experiment.setLog2PseudoSPCThreshold(mean2 + timesSigmaForSpecificity * sigma2, replicate);
		experiment.setGaussianFit(gaussian, replicate);
	}

	private void calculateAndFitNoise(TurboIDExperiment experiment, File outputFileGaussian, double sigmaFactorToFit)
			throws IOException {

		final List<WeightedObservedPoint> points = getHistogramOfLog2PseudoSPC(experiment);
		if (points.isEmpty()) {
			return;
		}
		experiment.setExperimentalHistogram(points);
		final double[] fit = GaussianCurveFitter.create().fit(points);
		final double normalization = fit[0];
		final double mean = fit[1];
		final double sigma = fit[2];
		System.out.println(normalization + "\t" + mean + "\t" + sigma);

		// log.info("Now, I mirror the right part of the distribution from mean
		// to max and then, I do fitting again");
		// final List<WeightedObservedPoint> newPoints =
		// getNewPointsToFitMirroringAfterInflexionPoint(points);
		log.info("Now, I take only up to " + sigmaFactorToFit + "x sigma to fit again the data");
		final List<WeightedObservedPoint> newPoints = getNewPointsToFitUpTo1Sigma(points, mean,
				sigmaFactorToFit * sigma);

		final double[] fit2 = GaussianCurveFitter.create().fit(newPoints);
		final double normalization2 = fit2[0];
		final double mean2 = fit2[1];
		final double sigma2 = fit2[2];
		System.out.println(normalization2 + "\t" + mean2 + "\t" + sigma2);

		final Gaussian gaussian = new Gaussian(normalization2, mean2, sigma2);
		final FileWriter fw = new FileWriter(outputFileGaussian);
		fw.write("Gaussian fitting:\n");
		fw.write("Normalization=" + normalization2 + "\n");
		fw.write("Mean=" + mean2 + "\n");
		fw.write("Sigma=" + sigma2 + "\n\n");
		fw.write("log2 pseudoSPC" + "\t" + "observed" + "\t" + "used for fitting" + "\t" + "fitted" + "\n");

		for (int index = 0; index < points.size(); index++) {
			final WeightedObservedPoint point = points.get(index);

			final double x = point.getX();
			final double observed = point.getY();
			String observedMirroredString = "";
			if (newPoints.size() > index) {
				final WeightedObservedPoint point2 = newPoints.get(index);
				observedMirroredString = String.valueOf(point2.getY());
			}
			final double fitted = gaussian.value(x);
			fw.write(x + "\t" + observed + "\t" + observedMirroredString + "\t" + fitted + "\n");
		}
		fw.close();
		log.info("File written at: " + outputFileGaussian.getAbsolutePath());

		experiment.setLog2PseudoSPCThreshold(mean2 + timesSigmaForSpecificity * sigma2);
		experiment.setGaussianFit(gaussian);
	}

	private List<WeightedObservedPoint> getNewPointsToFitUpTo1Sigma(List<WeightedObservedPoint> points, double mean,
			double sigma) {
		final List<WeightedObservedPoint> newPoints = new ArrayList<WeightedObservedPoint>();

		final double thresholdToFit = mean + 1 * sigma;
		for (int i = 0; i < points.size(); i++) {
			if (points.get(i).getX() <= thresholdToFit) {

				newPoints.add(i, points.get(i));
			}
		}
		return newPoints;
	}

	private List<WeightedObservedPoint> getHistogramOfLog2PseudoSPC(TurboIDExperiment experiment, Replicate replicate) {
		final TDoubleArrayList data = new TDoubleArrayList();
		final Collection<ProteinFromTurboID> proteins = experiment.values();
		for (final ProteinFromTurboID protein : proteins) {
			final double[] pseudoSPCs = protein.getPseudoSpecCountsWithNormalizedIntensities(replicate).values();
			for (final double pseudoSPC : pseudoSPCs) {
				if (!Double.isNaN(pseudoSPC) && Double.compare(pseudoSPC, 0.0) != 0) {
					final double log2 = Maths.log(pseudoSPC, 2);
					data.add(log2);
				}
			}
		}
		if (data.isEmpty()) {
			return Collections.emptyList();
		}
		final int numBins = 35;
		final int[] histogram = Histogram.calcHistogram(data, numBins);
		final double binSize = Histogram.getBinSize(data, numBins);
		final List<WeightedObservedPoint> ret = new ArrayList<WeightedObservedPoint>();

		double point = data.min();
		int index = 0;
		while (point < data.max()) {
			if (index == histogram.length) {
				break;
			}
			final int frecuency = histogram[index++];
			ret.add(new WeightedObservedPoint(1, point, frecuency));
			point += binSize;
		}
		return ret;
	}

	private List<WeightedObservedPoint> getHistogramOfLog2PseudoSPC(TurboIDExperiment experiment) {
		final TDoubleArrayList data = new TDoubleArrayList();
		final Collection<ProteinFromTurboID> proteins = experiment.values();
		for (final ProteinFromTurboID protein : proteins) {
			for (final Replicate replicate : Replicate.values(experiment.getFraction())) {
				final double[] pseudoSPCs = protein.getPseudoSpecCountsWithNormalizedIntensities(replicate).values();
				for (final double pseudoSPC : pseudoSPCs) {
					if (!Double.isNaN(pseudoSPC) && Double.compare(pseudoSPC, 0.0) != 0) {
						final double log2 = Maths.log(pseudoSPC, 2);
						data.add(log2);
					}
				}
			}
		}
		if (data.isEmpty()) {
			return Collections.emptyList();
		}
		final int numBins = 35;
		final int[] histogram = Histogram.calcHistogram(data, numBins);
		final double binSize = Histogram.getBinSize(data, numBins);
		final List<WeightedObservedPoint> ret = new ArrayList<WeightedObservedPoint>();

		double point = data.min();
		int index = 0;
		while (point < data.max()) {
			if (index == histogram.length) {
				break;
			}
			final int frecuency = histogram[index++];
			ret.add(new WeightedObservedPoint(1, point, frecuency));
			point += binSize;
		}
		return ret;
	}

	private TurboIDExperiment readTurboIDExperiment(File file, String[] sheets, TurboIDFraction fraction)
			throws IOException {
		final Map<String, String> genesByACC = new THashMap<String, String>();
		final TurboIDExperiment turboIDExperiment = new TurboIDExperiment(fraction);
		final ExcelReader reader = new ExcelReader(file, 0, 6);
		for (final String sheetName : sheets) {
			final Replicate replicate = replicateBySheets.get(sheetName);
			final Set<String> genes = new THashSet<String>();
			final int sheetIndex = reader.getWorkbook().getSheetIndex(sheetName);
			if (sheetIndex < 0) {
				throw new IllegalArgumentException(
						sheetName + " sheet is not found in input file " + file.getAbsolutePath());
			}
			final int columnIndexForACC = reader.getColumnIndex(sheetIndex, ACCESSION);
			final int columnIndexForDescription = reader.getColumnIndex(sheetIndex, DESCRIPTION);
			final int columnIndexForSPC = reader.getColumnIndex(sheetIndex, SPC);
			int numRow = 1;
			while (true) {
				try {
					final String acc = reader.getStringValue(sheetIndex, numRow, columnIndexForACC);
					if (acc == null) {
						break;
					}
					if (FastaParser.isReverse(acc) || FastaParser.isContaminant(acc)) {
						continue;
					}
					final String description = reader.getStringValue(sheetIndex, numRow, columnIndexForDescription);
					int spc = 0;
					if (columnIndexForSPC >= 0) {
						spc = Double.valueOf(reader.getStringValue(sheetIndex, numRow, columnIndexForSPC)).intValue();
					}
					String gene = parseForExcelConflictingGenes(FastaParser.getGeneFromFastaHeader(description));
					boolean problematic = false;
					while (genes.contains(gene)) {
						if (genesByACC.containsKey(acc)) {
							gene = genesByACC.get(acc);
							break;
						}
						gene = gene + "-" + acc;
						problematic = true;
					}
					if (problematic) {
						genesByACC.put(acc, gene);
					}
					genes.add(gene);
					boolean isTransmembrane = false;
					final int columnIndex = reader.getColumnIndex(sheetIndex, TRANSMEMBRANE);
					if (columnIndex >= 0) {
						final String transmembraneString = reader.getStringValue(sheetIndex, numRow, columnIndex);
						if (transmembraneString != null && !"".equals(transmembraneString)) {
							isTransmembrane = Boolean.valueOf(transmembraneString);
						}
					}
					ProteinFromTurboID protein = null;
					if (turboIDExperiment.containsKey(acc)) {
						protein = turboIDExperiment.get(acc);
					} else {
						protein = new ProteinFromTurboID(acc, gene, isTransmembrane);
						turboIDExperiment.put(protein.getAcc(), protein);
					}
					protein.setSpc(spc, replicate);
				} finally {
					numRow++;
				}
			}
		}
		log.info(turboIDExperiment.size() + " proteins in total");
		// then, read all the intensities

		for (final String sheetName : sheets) {
			final int sheetIndex = reader.getWorkbook().getSheetIndex(sheetName);
			if (sheetIndex < 0) {
				throw new IllegalArgumentException(
						"Sheet " + sheetName + " not found in excel file " + file.getAbsolutePath());
			}
			final int columnIndexForACC = reader.getColumnIndex(sheetIndex, ACCESSION);

			int numRow = 2;
			while (true) {
				try {
					final String acc = reader.getStringValue(sheetIndex, numRow, columnIndexForACC);
					if (acc == null) {
						break;
					}
					if (!turboIDExperiment.containsKey(acc)) {
						continue;
					}
					final ProteinFromTurboID protein = turboIDExperiment.get(acc);

					///////////////////////////////
					final Set<TurboID_Channel_Ori> notFound2 = new THashSet<TurboID_Channel_Ori>();
					int counter = 0;
					for (final TurboID_Channel_Ori replicateID : TurboID_Channel_Ori.values()) {
						final int columnIndex = reader.getColumnIndex(sheetIndex, replicateID.name());
						if (columnIndex == -1) {
							notFound2.add(replicateID);
							continue;
						}
						counter++;
						final String stringValue = reader.getNumberValue(sheetIndex, numRow, columnIndex);
						Double intensity = Double.NaN;
						if (stringValue != null && !stringValue.equals("NA")) {
							try {
								intensity = Double.valueOf(stringValue);
							} catch (final NumberFormatException e) {
								throw new IllegalArgumentException("Error reading on sheet " + sheetName + " line "
										+ (numRow + 1) + " column " + (columnIndex + 1));
							}
						}
						final TurboIDExperimentType turboIDExperimenType = replicateID.getExpType();
						protein.addOriginalIntensity(intensity, replicateID, turboIDExperimenType);
					}
					if (counter != TurboID_Channel_Ori.values().length / 6) {
						for (final TurboID_Channel_Ori string : notFound2) {
							log.info(string);
						}
						throw new IllegalArgumentException("Some columns were not found");
					}
					/////////////////////////////
					if (!normalizeChannelsManually) {
						final Set<TurboID_Channel_Norm> notFound = new THashSet<TurboID_Channel_Norm>();
						counter = 0;
						final int index = 0;
						for (final TurboID_Channel_Norm channel : TurboID_Channel_Norm.values()) {

							final int columnIndex = reader.getColumnIndex(sheetIndex, channel.name());
							if (columnIndex == -1) {
								notFound.add(channel);
								continue;
							}
							counter++;
							final String stringValue = reader.getNumberValue(sheetIndex, numRow, columnIndex);
							Double intensity = Double.NaN;
							if (stringValue != null && !stringValue.equals("NA")) {
								intensity = Double.valueOf(stringValue);
							}
							final TurboIDExperimentType turboIDExperimenType = channel.getExpType();
							protein.addNormalizedIntensity(intensity, channel, turboIDExperimenType);
						}

						if (counter != TurboID_Channel_Norm.values().length / 2) {
							for (final TurboID_Channel_Norm string : notFound) {
								log.info(string);
							}
							throw new IllegalArgumentException("Some columns were not found");
						}
					}
				} finally {
					numRow++;
				}
			}

		}
		if (minimumIntensityThreshold != null) {
			int numDiscarded = 0;
			final Set<String> toRemove = new THashSet<String>();
			for (final String key : turboIDExperiment.keySet()) {
				final ProteinFromTurboID protein = turboIDExperiment.get(key);
				boolean valid = false;
				for (final TurboID_Channel_Ori channelOri : TurboID_Channel_Ori.values()) {
					final Double intensity = protein.getOriginalIntensities(channelOri.getExpType()).get(channelOri);
					if (intensity >= minimumIntensityThreshold) {
						valid = true;
						break;
					}
				}
				if (!valid) {
					// remove protein
					numDiscarded++;
					toRemove.add(key);
				}
			}
			for (final String key : toRemove) {
				turboIDExperiment.remove(key);
			}
			log.info(numDiscarded
					+ " proteins were discarded because they didn't have any channel intensity above threshold "
					+ minimumIntensityThreshold);
		}
		if (normalizeChannelsManually) {
			for (final ProteinFromTurboID protein : turboIDExperiment.values()) {

				int index = 0;
				for (final TurboID_Channel_Norm channel : TurboID_Channel_Norm.values()) {
					final TurboID_Channel_Ori channelOri = TurboID_Channel_Ori.values()[index++];
					final double sumChannelIntensity = turboIDExperiment.getSumIntensities(channelOri);
					final Double intensity = protein.getOriginalIntensities(channel.getExpType()).get(channelOri);
					final double normalizedIntensity = intensity / sumChannelIntensity;
					final TurboIDExperimentType turboIDExperimenType = channel.getExpType();
					protein.addNormalizedIntensity(normalizedIntensity, channel, turboIDExperimenType);
				}
			}
		}
		return turboIDExperiment;

	}

	protected static boolean isInExclusionList(String acc) {
		if ("P22629".equalsIgnoreCase(acc)) {
			return true;
		}
		if (acc.toLowerCase().contains("v5tag")) {
			return true;
		}
		return false;
	}

	private void loadUniprotAnnotations(Collection<String> accs) {
		log.info("Loading Uniprot annotations for " + accs.size() + " proteins");
		final UniprotProteinRetriever upr = new UniprotProteinRetriever(null, uniprotReleasesFolder, true);
		upr.getAnnotatedProteins(accs);
		log.info("annotations loaded");
	}

	private String getOutputFileName(String prefix, boolean onlyTM, boolean onlyNonTM, boolean onlyComplete,
			boolean ignoreMissingValues, boolean useOfNormalizedIntensities, boolean normalizeToAverage,
			boolean normalizeToTurboIDChannelAverage, boolean log2Applied) {
		final String v5 = "";
		// if (normalizingByV5Tag) {
		// v5 = "_V5";
		// }
		final String mixChannel = "";
		// if (normalizingByMixChannel) {
		// mixChannel = "_mix";
		// }
		String tmString = "_all";
		if (onlyTM || onlyNonTM || onlyComplete) {
			if (onlyTM) {
				tmString = "_TM";
			}
			if (onlyNonTM) {
				tmString = "_nonTM";
			}
			if (onlyComplete) {
				tmString = "_complete";
			}
		}

		String ignoreMissingValuesString = "";
		if (!ignoreMissingValues) {
			ignoreMissingValuesString = "_withNA";
		}
		String useNormalizedIntensitiesString = "_norm";
		if (!useOfNormalizedIntensities) {
			useNormalizedIntensitiesString = "_ori";
		}
		String normalizeToAverageString = "";
		if (normalizeToAverage) {
			normalizeToAverageString = "_avgNorm";
		}
		String normalizeToTurboIDChannelAverageString = "";
		if (normalizeToTurboIDChannelAverage) {
			normalizeToTurboIDChannelAverageString = "_turboIDNorm";
		}
		String log2AppliedString = "";
		if (log2Applied) {
			log2AppliedString = "_log2";
		}
		return prefix + tmString + useNormalizedIntensitiesString + mixChannel + v5 + normalizeToAverageString
				+ normalizeToTurboIDChannelAverageString + log2AppliedString + ignoreMissingValuesString + ".txt";
	}

	private void setupChannels() {
		replicateBySheets.put(sheetsTMT11Nu[0], Replicate.A);
		if (sheetsTMT11Nu.length > 1) {
			replicateBySheets.put(sheetsTMT11Nu[1], Replicate.B);
		}
		if (sheetsTMT11Nu.length > 2) {
			replicateBySheets.put(sheetsTMT11Nu[2], Replicate.Brerun);
		}
		if (sheetsTMT11Nu.length > 3) {
			replicateBySheets.put(sheetsTMT11Nu[3], Replicate.Arerun);
		}

		replicateBySheets.put(sheetsTMT11Cy[0], Replicate.CyA);
		if (sheetsTMT11Cy.length > 1) {
			replicateBySheets.put(sheetsTMT11Cy[1], Replicate.CyB);
		}
	}

	private String parseForExcelConflictingGenes(String gene) {
		if (gene == null) {
			gene = "N/A";
		}
		if (gene.toLowerCase().startsWith("march")) {
			gene = "'" + gene;
		}
		if (gene.toLowerCase().startsWith("marc")) {
			gene = "'" + gene;
		}
		if (gene.toLowerCase().startsWith("septin")) {
			gene = "'" + gene;
		}
		if (gene.toLowerCase().startsWith("sept")) {
			gene = "'" + gene;
		}
		return gene;
	}

	private void exportSpecificity(TurboIDExperiment turboIDExperiment, File outputFile) throws IOException {
		final FileWriter fw = new FileWriter(outputFile);
		final List<ProteinFromTurboID> proteins = turboIDExperiment
				.getProteinsSorted(TurboIDExperiment.getComparatorByTotalSPC(false));
		fw.write("accession" + "\t" + "gene" + "\t" + "transmembrane" + "\t" + "total SPCs" + "\t");
		// spc per replicate
		for (final Replicate replicate : Replicate.values(turboIDExperiment.getFraction())) {
			fw.write("SPC_" + replicate.name() + "\t");
		}
		for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaits()) {

			fw.write("FDR_" + bait.name() + "\t");

		}
		// for (final TurboIDExperimentType bait :
		// TurboIDExperimentType.getBaits()) {
		// fw.write(bait.name() + "\t");
		// }
		// for (final Replicate replicate : Replicate.values(fraction)) {
		// for (final TurboIDExperimentType bait :
		// TurboIDExperimentType.getBaits()) {
		// fw.write(bait.name() + "_" + replicate.name() + "\t");
		// }
		// fw.write("ALL_" + replicate.name() + "\t");
		// }

		for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaits()) {
			for (final Replicate replicate : Replicate.values(turboIDExperiment.getFraction())) {
				fw.write("AVG_pseudoSPC_" + bait.name() + "_" + replicate.name() + "\t");
			}
		}
		// for (final TurboIDExperimentType bait :
		// TurboIDExperimentType.getBaits()) {
		// final TurboID_Channel_Norm[] channels =
		// TurboID_Channel_Norm.values();
		// for (final TurboID_Channel_Norm channel : channels) {
		// if (channel.getExpType() == bait) {
		// fw.write("SPC_" + channel.name() + "\t");
		// }
		// }
		//
		// }
		// intensities
		for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaits()) {
			for (final Replicate replicate : Replicate.values(turboIDExperiment.getFraction())) {
				fw.write("AVG_normInt_" + bait.name() + "_" + replicate.name() + "\t");
			}
		}
		fw.write("\n");

		for (final ProteinFromTurboID protein : proteins) {
			if (protein.getAcc().equals("Q9JM61")) {
				log.info(protein);
			}
			fw.write(protein.getAcc() + "\t" + protein.getGene() + "\t" + protein.isTransmembrane() + "\t"
					+ protein.getSumSPCAcrossReplicates() + "\t");
			// spc per replicate
			for (final Replicate replicate : Replicate.values(turboIDExperiment.getFraction())) {
				fw.write(protein.getSpc(replicate) + "\t");
			}
			// fdr per replicate per bait
			for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaits()) {
				final Gaussian fit = turboIDExperiment.getGaussianFit();
				final List<WeightedObservedPoint> experimentalHistogram = turboIDExperiment.getExperimentalHistogram();
				fw.write(protein.getLocalFDR(bait, fit, experimentalHistogram) + "\t");
			}
			// for (final TurboIDExperimentType bait :
			// TurboIDExperimentType.getBaits()) {
			// fw.write(protein.isSpecificInBothExperiments(bait,
			// turboIDExperiment.getLog2PseudoSPCThreshold(Replicate.A),
			// turboIDExperiment.getLog2PseudoSPCThreshold(Replicate.B)) +
			// "\t");
			// }

			// for (final Replicate replicate : Replicate.values(fraction)) {
			// boolean allInRep = true;
			// for (final TurboIDExperimentType bait :
			// TurboIDExperimentType.getBaits()) {
			// final boolean specific = protein.isSpecific(bait, replicate,
			// turboIDExperiment.getLog2PseudoSPCThreshold(replicate));
			// allInRep = allInRep && specific;
			// fw.write(specific + "\t");
			// }
			// fw.write(allInRep + "\t");
			// }

			for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaits()) {
				for (final Replicate replicate : Replicate.values(turboIDExperiment.getFraction())) {
					final TObjectDoubleHashMap<TurboID_Channel_Norm> pseudoSpecCountsWithNormalizedIntensities = protein
							.getPseudoSpecCountsWithNormalizedIntensities(replicate);
					final TDoubleArrayList toAverage = new TDoubleArrayList();
					for (final TurboID_Channel_Norm channel : pseudoSpecCountsWithNormalizedIntensities.keySet()) {
						if (channel.getExpType() == bait) {
							toAverage.add(pseudoSpecCountsWithNormalizedIntensities.get(channel));
						}
					}
					final double mean = Maths.mean(toAverage);
					if (Double.isNaN(mean)) {
						fw.write(mean + "\t");
					} else {
						fw.write(formatter.format(mean) + "\t");
					}
				}
			}
			// for (final TurboIDExperimentType expType :
			// TurboIDExperimentType.getBaits()) {
			// final TurboID_Channel_Norm[] channels =
			// TurboID_Channel_Norm.values();
			// for (final TurboID_Channel_Norm channel : channels) {
			// if (channel.getExpType() == expType) {
			// final TObjectDoubleHashMap<TurboID_Channel_Norm>
			// pseudoSpecCountsWithNormalizedIntensities = protein
			// .getPseudoSpecCountsWithNormalizedIntensities(channel.getReplicate());
			// fw.write(formatter.format(pseudoSpecCountsWithNormalizedIntensities.get(channel))
			// + "\t");
			// }
			// }
			// }
			// intensities
			for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaits()) {
				final TObjectDoubleHashMap<TurboID_Channel_Norm> normalizedIntensitiesPerChannel = protein
						.getNormalizedIntensities(bait);
				for (final Replicate replicate : Replicate.values(turboIDExperiment.getFraction())) {
					final TDoubleArrayList toAverage = new TDoubleArrayList();
					for (final TurboID_Channel_Norm channel : normalizedIntensitiesPerChannel.keySet()) {
						if (channel.getReplicate() == replicate) {
							toAverage.add(normalizedIntensitiesPerChannel.get(channel));
						}
					}
					final double mean = Maths.mean(toAverage);
					fw.write(mean + "\t");
				}
			}
			fw.write("\n");
		}
		fw.close();
		log.info("File written at: " + outputFile.getAbsolutePath());
	}
}
