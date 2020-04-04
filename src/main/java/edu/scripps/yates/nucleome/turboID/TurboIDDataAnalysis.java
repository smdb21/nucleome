package edu.scripps.yates.nucleome.turboID;

import static org.junit.Assert.fail;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
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
import smile.netlib.NLMatrix;

public class TurboIDDataAnalysis {
	private final static Logger log = Logger.getLogger(TurboIDDataAnalysis.class);
	private static final String ACCESSION = "accession";
	private static final String TRANSMEMBRANE = "Transmembrane";
	private static final String DESCRIPTION = "description";
	private static final String SPC = "spec count";
	private static final String PSM = "psm";
	private static final int GO_TERM_DISTANCE = 3;
	private final File excelFile_TMT11_TurboID_Nu = new File(
			"Z:\\share\\Salva\\data\\4D_Nucleome\\TurboID\\input\\NE4_Nu_T11aT11b_QtS03_TagM_QtRaw.xlsx");
	private final String[] sheetsTMT11Nu = { "T11a_QtS03_TagM_trimNL", "T11b_QtS03_TagM_trimNL", "1114LT_Nu_Tb_rerun" };
	private final File excelFile_TMT11_TurboID_Cy = new File(
			"Z:\\share\\Salva\\data\\4D_Nucleome\\TurboID\\input\\RawData_Cytoplam_TMT_Run1andRun2.xlsx");
	private final String[] sheetsTMT11Cy = { "Cytoplasmic_Run1", "Cytoplasmic_Run2" };
	private final File excelFile_TMT8_EMD_IP = new File(
			"Z:\\share\\Salva\\data\\4D_Nucleome\\TMT8_EMD_IP\\4DN_0518LT_IP_T8_Qt_pro_Venn.xlsx");
	private final String sheet_TMT8_EMD_IP = "raw_0518LT_IP_T8_Qt_pro_csv";

	private final File outputFolder = excelFile_TMT11_TurboID_Nu.getParentFile();
	protected static final File uniprotReleasesFolder = new File("Z:\\share\\Salva\\data\\uniprotKB");
	private final AnnotationsUtil annotationsUtil = new AnnotationsUtil(uniprotReleasesFolder);
	private final boolean[] yesno = { true, false };
	// private boolean normalizedByMixChannel;
	// private boolean normalizeByV5Tag;
	private final THashMap<String, Replicate> replicateBySheets = new THashMap<String, Replicate>();
	private boolean normalizeChannelsManually;
	private double timesSigmaForSpecificity;
	private final DecimalFormat formatter = new DecimalFormat("#.#");

	// *************************
	// FRACTION TO ANALYZE
	// *************************
	private final boolean analyzeNucleusFraction = true;
	private final boolean analyzeCytoplasmFraction = false;
	private TurboIDFraction fraction;
	// *************************

	public static void main(String[] args) {
		new TurboIDDataAnalysis().run();
	}

	public void run() {
		try {
			// take first the TMT8 EMD IP experiment
			final IPExperiment ipExperiment = createIPExperiment();

			setupChannels();
			final boolean applyLog2 = true;

			// for (final boolean useNormalizedIntensities : yesno) {
			final boolean onlyTM = false;
			final boolean onlyNonTM = false;
			final boolean useNormalizedIntensities = false;
			final boolean ignoreMissingValues = true;
			// normalizations:
			final boolean normalizedByMixChannel = false; // discarded, always
															// FALSE
			final boolean normalizeByV5Tag = false; // discarded, always FALSE
			final boolean normalizeToAverage = false;
			final boolean normalizeToTurboIDChannelAverage = false;
			final boolean onlySignificantlySpecific = true;
			normalizeChannelsManually = true; // normalize by the sum of the
												// channel
			timesSigmaForSpecificity = 2;

			TurboIDExperiment turboIDExperiment = null;
			String fileoutputPrefix = null;
			int minPerGroup = 0;
			if (analyzeNucleusFraction) {
				turboIDExperiment = readTurboIDExperiment(excelFile_TMT11_TurboID_Nu, sheetsTMT11Nu);
				fileoutputPrefix = "Nu";
				fraction = TurboIDFraction.NU;
				minPerGroup = 4;
			} else if (analyzeCytoplasmFraction) {
				turboIDExperiment = readTurboIDExperiment(excelFile_TMT11_TurboID_Cy, sheetsTMT11Cy);
				fileoutputPrefix = "Cy";
				fraction = TurboIDFraction.CY;
				minPerGroup = 4;
			}
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
			final File outputFile = new File(outputFolder.getAbsolutePath() + File.separator
					+ getOutputFileName(fileoutputPrefix, onlyTM, onlyNonTM, ignoreMissingValues,
							useNormalizedIntensities, normalizeToAverage, normalizeToTurboIDChannelAverage, applyLog2));

			final double sigmaFactor = 0.25;
			for (final Replicate replicate : Replicate.values(fraction)) {
				final File fileGaussianRep = new File(outputFile.getParentFile().getAbsoluteFile() + File.separator
						+ FilenameUtils.getBaseName(outputFile.getAbsolutePath()) + "_GaussFit_" + replicate.name()
						+ ".txt");
				calculateAndFitNoise(turboIDExperiment, replicate, fileGaussianRep, sigmaFactor);
			}
			final File fileGaussian = new File(outputFile.getParentFile().getAbsoluteFile() + File.separator
					+ FilenameUtils.getBaseName(outputFile.getAbsolutePath()) + "_GaussFit.txt");
			calculateAndFitNoise(turboIDExperiment, fileGaussian, sigmaFactor);

			turboIDExperiment.exportToFileForClustering(outputFile, onlyTM, onlyNonTM, ignoreMissingValues,
					useNormalizedIntensities);
			if (analyzeNucleusFraction) {
				final File outputFile2 = new File(outputFile.getParentFile().getAbsoluteFile() + File.separator
						+ FilenameUtils.getBaseName(outputFile.getAbsolutePath()) + "_vs_TMT8Nu.txt");
				exportWholeNucleiFractionFromIP_vs_TurboIDOnlyFromTurboID(turboIDExperiment, ipExperiment, outputFile2,
						applyLog2, useNormalizedIntensities, onlyTM, onlyNonTM, onlySignificantlySpecific);
				// //
				//
				// final File outputFile3 = new
				// File(outputFile.getParentFile().getAbsoluteFile() +
				// File.separator
				// + FilenameUtils.getBaseName(outputFile.getAbsolutePath()) +
				// "_vs_TMT8Nu_correl_matrix.txt");
				// exportWholeNucleiFractionFromIP_vs_TurboIDOnlyFromTurboID_correlationMatrix(turboIDExperiment,
				// ipExperiment, outputFile3, applyLog2,
				// useNormalizedIntensities, onlyTM, onlyNonTM,
				// onlySignificantlySpecific);

				//
				final File outputFile4 = new File(outputFile.getParentFile().getAbsoluteFile() + File.separator
						+ FilenameUtils.getBaseName(outputFile.getAbsolutePath()) + "_error_distributions.txt");
				exportReplicateErrorDistribution(turboIDExperiment, outputFile4, applyLog2, useNormalizedIntensities);
				//
				final File outputFile5 = new File(outputFile.getParentFile().getAbsoluteFile() + File.separator
						+ FilenameUtils.getBaseName(outputFile.getAbsolutePath()) + "_specific.txt");
				exportSpecificity(turboIDExperiment, outputFile5);
			}
			// calculate ANOVA per protein comparing the baits

			final double anovaPValueThreshold = 0.05;

			// now the same but with distribSPC
			final File outputFile8 = new File(outputFile.getParentFile().getAbsoluteFile() + File.separator
					+ fileoutputPrefix + "_anova_log_distribSPC.txt");
			calculateAnovaPerProteinComparingBaits(turboIDExperiment, outputFile8, true, true, false,
					anovaPValueThreshold, minPerGroup);
			// now the same but with distribSPC
			final File outputFile9 = new File(outputFile.getParentFile().getAbsoluteFile() + File.separator
					+ fileoutputPrefix + "_anova_log_distribInt.txt");
			calculateAnovaPerProteinComparingBaits(turboIDExperiment, outputFile9, true, false, true,
					anovaPValueThreshold, minPerGroup);
			final File outputFile10 = new File(outputFile.getParentFile().getAbsoluteFile() + File.separator
					+ fileoutputPrefix + "_anova_distribInt.txt");
			calculateAnovaPerProteinComparingBaits(turboIDExperiment, outputFile10, false, false, true,
					anovaPValueThreshold, minPerGroup);
		} catch (final Exception e) {
			e.printStackTrace();
			fail();
		}

	}

	private void calculateAnovaPerProteinComparingBaits(TurboIDExperiment turboIDExperiment, File outputFile,
			boolean applyLog, boolean distribSPC, boolean distribIntensity, double anovaPValueThreshold,
			int minPerGroup) throws IOException {
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
		fw.write(prefix + " ANOVA p-value " + "\t");
		for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaits()) {
			fw.write(prefix + " " + bait.name() + "\t");
		}
		for (final TurboIDExperimentType bait : TurboIDExperimentType.values()) {
			for (final Replicate replicate : Replicate.values(fraction)) {
				for (final TurboID_Channel_Norm channel : TurboID_Channel_Norm.values()) {
					if (channel.getExpType() == bait && channel.getReplicate() == replicate) {

						fw.write(prefix + " " + bait.name() + "_" + channel.name() + "\t");

					}
				}
			}
		}
		fw.write("\n");
		for (final ProteinFromTurboID protein : proteins) {
			if (protein.getAcc().equals("Q9JM61")) {
				log.info(protein);
			}
			fw.write(protein.getAcc() + "\t" + protein.getGene() + "\t" + annotationsUtil.getDescription(protein) + "\t"
					+ annotationsUtil.getTransmembraneRegion(protein) + "\t" + annotationsUtil.isNucleus(protein) + "\t"
					+ annotationsUtil.isDNABinding(protein) + "\t" + annotationsUtil.isTranscriptionFactor(protein)
					+ "\t" + annotationsUtil.isRNABinding(protein) + "\t" + annotationsUtil.isHeterochromatin(protein)
					+ "\t" + protein.getSumSPCAcrossReplicates() + "\t");
			// spc per replicate
			for (final Replicate replicate : Replicate.values(fraction)) {
				fw.write(protein.getSpc(replicate) + "\t");
			}
			// column saying whether contains enough data to do the pvalue or
			// not

			final boolean isvalid = protein.isAnovaValid(distribSPC, distribIntensity, fraction, minPerGroup);
			fw.write(isvalid + "\t");
			final double anovaPValue = protein.getAnovaPValueOverBaits(applyLog, distribSPC, distribIntensity,
					fraction);
			fw.write(anovaPValue + "\t");
			final TObjectDoubleHashMap<TurboID_Channel_Norm> distributedIntensitiesWithNormalizedIntensities = protein
					.getDistributedIntensitiesWithNormalizedIntensities();
			for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaits()) {
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
						value = Double.NaN;
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

			for (final TurboIDExperimentType bait : TurboIDExperimentType.values()) {

				for (final Replicate replicate : Replicate.values(fraction)) {
					for (final TurboID_Channel_Norm channel : TurboID_Channel_Norm.values()) {
						if (channel.getExpType() == bait && channel.getReplicate() == replicate) {
							double value = Double.NaN;
							if (distribSPC) {
								final TObjectDoubleHashMap<TurboID_Channel_Norm> pseudoSpecCountsWithNormalizedIntensities = protein
										.getPseudoSpecCountsWithNormalizedIntensities(replicate);
								if (pseudoSpecCountsWithNormalizedIntensities.containsKey(channel)) {
									value = pseudoSpecCountsWithNormalizedIntensities.get(channel);
								}
							} else if (distribIntensity) {

								if (distributedIntensitiesWithNormalizedIntensities.containsKey(channel)) {
									value = distributedIntensitiesWithNormalizedIntensities.get(channel);
								}
							} else {
								value = protein.getNormalizedIntensities(bait).get(channel);
							}

							if (applyLog && !Double.isNaN(value)) {
								final double log2 = Maths.log(value, 2);
								if (Double.isInfinite(log2)) {
									log.info("asdf");
								}
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

	private void exportSpecificity(TurboIDExperiment turboIDExperiment, File outputFile) throws IOException {
		final FileWriter fw = new FileWriter(outputFile);
		final List<ProteinFromTurboID> proteins = turboIDExperiment
				.getProteinsSorted(TurboIDExperiment.getComparatorByTotalSPC(false));
		fw.write("accession" + "\t" + "gene" + "\t" + "transmembrane" + "\t" + "total SPCs" + "\t");
		// spc per replicate
		for (final Replicate replicate : Replicate.values(fraction)) {
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
			for (final Replicate replicate : Replicate.values(fraction)) {
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
			for (final Replicate replicate : Replicate.values(fraction)) {
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
			for (final Replicate replicate : Replicate.values(fraction)) {
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
				for (final Replicate replicate : Replicate.values(fraction)) {
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
				for (final Replicate replicate : Replicate.values(fraction)) {
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

	private void exportSpecificityOverlapsBetweenReps(TurboIDExperiment turboIDExperiment, File outputFile)
			throws IOException {
		final FileWriter fw = new FileWriter(outputFile);
		final List<ProteinFromTurboID> proteins = turboIDExperiment
				.getProteinsSorted(TurboIDExperiment.getComparatorByTotalSPC(false));
		fw.write("accession" + "\t" + "gene" + "\t" + "transmembrane" + "\t" + "total SPCs" + "\t");
		// spc per replicate
		for (final Replicate replicate : Replicate.values(fraction)) {
			fw.write("SPC_" + replicate.name() + "\t");
		}

		fw.close();
		log.info("File written at: " + outputFile.getAbsolutePath());
	}

	private void calculateAndFitNoise(TurboIDExperiment experiment, Replicate replicate, File outputFileGaussian,
			double sigmaFactorToFit) throws IOException {

		final List<WeightedObservedPoint> points = getHistogramOfLog2PseudoSPC(experiment, replicate);
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

	private List<WeightedObservedPoint> getNewPointsToFitMirroringAfterInflexionPoint(
			List<WeightedObservedPoint> points) {
		final List<WeightedObservedPoint> newPoints = new ArrayList<WeightedObservedPoint>();

		final double max = getYMax(points);

		int j = 1;
		int inflexionPointIndex = -1;
		for (int i = 0; i < points.size(); i++) {
			if (inflexionPointIndex != -1
					|| (i > 0 && points.get(i).getY() > max * 0.8 && points.get(i).getY() < points.get(i - 1).getY())) {
				if (inflexionPointIndex == -1) {
					inflexionPointIndex = i;
				}
				j++;
				double newY = 0.0;
				if (j <= inflexionPointIndex) {
					newY = points.get(inflexionPointIndex - j).getY();
				}
				final double x = points.get(i).getX();
				newPoints.add(i, new WeightedObservedPoint(1.0, x, newY));
			} else {
				newPoints.add(points.get(i));
			}
		}
		return newPoints;
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

	private double getYMax(List<WeightedObservedPoint> points) {
		Double max = -Double.MAX_VALUE;
		for (final WeightedObservedPoint weightedObservedPoint : points) {
			if (max < weightedObservedPoint.getY()) {
				max = weightedObservedPoint.getY();
			}
		}
		return max;
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
			for (final Replicate replicate : Replicate.values(fraction)) {
				final double[] pseudoSPCs = protein.getPseudoSpecCountsWithNormalizedIntensities(replicate).values();
				for (final double pseudoSPC : pseudoSPCs) {
					if (!Double.isNaN(pseudoSPC) && Double.compare(pseudoSPC, 0.0) != 0) {
						final double log2 = Maths.log(pseudoSPC, 2);
						data.add(log2);
					}
				}
			}
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

	private void exportReplicateErrorDistribution(TurboIDExperiment turboIDExperiment, File outputFile,
			boolean applyLog2, boolean useNormalizedIntensities) throws IOException {
		final List<ProteinFromTurboID> proteins = turboIDExperiment
				.getProteinsSorted(TurboIDExperiment.getComparatorByTotalSPC(false));
		final FileWriter fw = new FileWriter(outputFile);
		fw.write("accession" + "\t" + "gene" + "\t" + "transmembrane" + "\t" + "total SPCs" + "\t");
		// spc per replicate
		for (final Replicate replicate : Replicate.values(fraction)) {
			fw.write("SPC_" + replicate.name() + "\t");
		}
		final TurboIDExperimentType[] expTypes = TurboIDExperimentType.values();
		for (final TurboIDExperimentType expType : expTypes) {
			if (expType.isBait()) {
				fw.write("AvgIntensity_" + expType.name() + "\t");
				fw.write("StdIntensity_" + expType.name() + "\t");
				fw.write("SEM_" + expType.name() + "\t");
				fw.write("Error_" + expType.name() + "\t");
				if (useNormalizedIntensities) {

					for (final TurboID_Channel_Norm channel : TurboID_Channel_Norm.values()) {
						if (channel.getExpType() == expType) {
							fw.write("Int_" + channel + "\t");
						}
					}
				} else {
					for (final TurboID_Channel_Ori channel : TurboID_Channel_Ori.values()) {
						if (channel.getExpType() == expType) {
							fw.write("Int_" + channel + "\t");
						}
					}
				}
				if (useNormalizedIntensities) {

					for (final TurboID_Channel_Norm channel : TurboID_Channel_Norm.values()) {
						if (channel.getExpType() == expType) {
							fw.write("pseudoSPC_" + channel + "\t");
						}
					}
				} else {
					for (final TurboID_Channel_Ori channel : TurboID_Channel_Ori.values()) {
						if (channel.getExpType() == expType) {
							fw.write("pseudoSPC_" + channel + "\t");
						}
					}
				}
			}
		}
		fw.write("\n");
		for (final ProteinFromTurboID protein : proteins) {
			fw.write(protein.getAcc() + "\t" + protein.getGene() + "\t" + protein.isTransmembrane() + "\t"
					+ protein.getSumSPCAcrossReplicates() + "\t");
			// spc per replicate
			for (final Replicate replicate : Replicate.values(fraction)) {
				fw.write(protein.getSpc(replicate) + "\t");
			}
			for (final TurboIDExperimentType expType : expTypes) {
				if (expType.isBait()) {

					final double avg = protein.getAvgOfIntensities(expType, useNormalizedIntensities);
					fw.write(avg + "\t");
					final double stdev = protein.getStandardDeviationOfIntensities(expType, useNormalizedIntensities);
					fw.write(stdev + "\t");
					final double error = protein.getStandardErrorOfMeasurementOfIntensities(expType,
							useNormalizedIntensities);
					fw.write(error + "\t");
					final double error2 = protein.getErrorMeasurement(expType, useNormalizedIntensities);
					fw.write(error2 + "\t");

					// intensities
					if (useNormalizedIntensities) {
						final TObjectDoubleHashMap<TurboID_Channel_Norm> normalizedIntensities = protein
								.getNormalizedIntensities(expType);
						for (final TurboID_Channel_Norm channel : TurboID_Channel_Norm.values()) {
							if (channel.getExpType() != expType) {
								continue;
							}
							if (normalizedIntensities.containsKey(channel)) {
								fw.write(normalizedIntensities.get(channel) + "");
							}
							fw.write("\t");

						}
					} else {
						final TObjectDoubleHashMap<TurboID_Channel_Ori> originalIntensities = protein
								.getOriginalIntensities(expType);
						for (final TurboID_Channel_Ori channel : TurboID_Channel_Ori.values()) {
							if (channel.getExpType() != expType) {
								continue;
							}
							if (originalIntensities.containsKey(channel)) {
								fw.write(originalIntensities.get(channel) + "");
							}
							fw.write("\t");

						}
					}
					// pseudoSPC
					if (useNormalizedIntensities) {
						for (final TurboID_Channel_Norm channel : TurboID_Channel_Norm.values()) {
							if (channel.getExpType() != expType) {
								continue;
							}
							final TObjectDoubleHashMap<TurboID_Channel_Norm> pseudoSpecCountsWithNormalizedIntensities = protein
									.getPseudoSpecCountsWithNormalizedIntensities(channel.getReplicate());
							fw.write("" + pseudoSpecCountsWithNormalizedIntensities.get(channel));
							fw.write("\t");
						}
					} else {
						for (final TurboID_Channel_Ori channel : TurboID_Channel_Ori.values()) {
							if (channel.getExpType() != expType) {
								continue;
							}

							final double d = protein.getPseudoSpecCountsWithOriginalIntensities(channel.getReplicate())
									.get(channel);
							fw.write("" + d);
							fw.write("\t");
						}
					}
				}
			}
			fw.write("\n");
		}

		fw.close();
		log.info("File written at: " + outputFile.getAbsolutePath());
	}

	private void exportWholeNucleiFractionFromIP_vs_TurboIDOnlyFromTurboID_correlationMatrix(
			TurboIDExperiment turboIDExperiment, IPExperiment ipExperiment, File outputFile, boolean applyLog2,
			boolean useNormalizedIntensities, boolean transmembraneOnly, boolean noTransmembraneOnly,
			boolean significantlySpecificOnly) throws IOException {

		final List<String> channels = new ArrayList<String>();
		channels.add("TMT8_Nu");
		if (useNormalizedIntensities) {
			for (final TurboID_Channel_Norm turboID_Channel_Norm : TurboID_Channel_Norm.values()) {
				if (turboID_Channel_Norm != TurboID_Channel_Norm.nmT1_mix_a
						&& turboID_Channel_Norm != TurboID_Channel_Norm.nmT1_mix_b) {
					channels.add(turboID_Channel_Norm.name());
				}
			}
		} else {
			for (final TurboID_Channel_Ori turboID_Channel_Ori : TurboID_Channel_Ori.values()) {
				if (turboID_Channel_Ori != TurboID_Channel_Ori.T1_mix_a
						&& turboID_Channel_Ori != TurboID_Channel_Ori.T1_mix_b) {
					channels.add(turboID_Channel_Ori.name());
				}
			}
		}
		final Collection<ProteinFromTurboID> proteins = turboIDExperiment.values();
		final NLMatrix matrix = new NLMatrix(channels.size(), channels.size(), Double.NaN);
		int row = 0;
		int col = 0;
		for (row = 0; row < channels.size(); row++) {
			for (col = row; col < channels.size(); col++) {
				final TDoubleArrayList array1 = new TDoubleArrayList();
				final TDoubleArrayList array2 = new TDoubleArrayList();
				for (final ProteinFromTurboID protein : proteins) {
					if (transmembraneOnly && !protein.isTransmembrane()) {
						continue;
					}
					if (noTransmembraneOnly && protein.isTransmembrane()) {
						continue;
					}

					Double val1 = null;
					Double val2 = null;
					if (useNormalizedIntensities) {
						final TObjectDoubleHashMap<TurboID_Channel_Norm> normalizedIntensitiesFromTurboID = new TObjectDoubleHashMap<TurboID_Channel_Norm>();
						for (final TurboIDExperimentType expType : TurboIDExperimentType.values()) {
							normalizedIntensitiesFromTurboID.putAll(protein.getNormalizedIntensities(expType));
						}
						// for array1, using row
						final TurboID_Channel_Norm channel1 = TurboID_Channel_Norm.getByName(channels.get(row));
						if (channel1 != null) {
							val1 = normalizedIntensitiesFromTurboID.get(channel1);
							if (applyLog2) {
								val1 = Maths.log(val1, 2);
							}

						} else {
							// if it is not returned by valueOf is because is
							// the TMT8_Nu
							final String acc = protein.getAcc();
							final ProteinFromIP proteinFromIP = ipExperiment.get(acc);
							if (proteinFromIP != null) {
								val1 = proteinFromIP.getNormalizedIntensity(IPExperimentType.NU);
								if (applyLog2) {
									val1 = Maths.log(val1, 2);
								}
							}

						}

						// for array2, using col
						final TurboID_Channel_Norm channel2 = TurboID_Channel_Norm.getByName(channels.get(col));
						if (channel2 != null) {
							val2 = normalizedIntensitiesFromTurboID.get(channel2);
							if (applyLog2) {
								val2 = Maths.log(val2, 2);
							}

						} else {
							// if it is not returned by valueOf is because is
							// the TMT8_Nu
							final String acc = protein.getAcc();
							final ProteinFromIP proteinFromIP = ipExperiment.get(acc);
							if (proteinFromIP != null) {
								val2 = proteinFromIP.getNormalizedIntensity(IPExperimentType.NU);
								if (applyLog2) {
									val2 = Maths.log(val2, 2);
								}
							}

						}
						// look if is significant in any bait in both
						// experiments
						if (significantlySpecificOnly) {
							if (channel1 == null || channel2 == null) {
								continue;
							}
							final boolean specificInBothExperiments = protein.isSpecificInBothExperiments(
									channel1.getExpType(), turboIDExperiment.getLog2PseudoSPCThreshold(Replicate.A),
									turboIDExperiment.getLog2PseudoSPCThreshold(Replicate.B));
							final boolean specificInBothExperiments2 = protein.isSpecificInBothExperiments(
									channel2.getExpType(), turboIDExperiment.getLog2PseudoSPCThreshold(Replicate.A),
									turboIDExperiment.getLog2PseudoSPCThreshold(Replicate.B));
							if (!specificInBothExperiments || !specificInBothExperiments2) {
								continue;
							}
						}

					} else {
						// using original intensities
						final TObjectDoubleHashMap<TurboID_Channel_Ori> originalIntensitiesFromTurboID = new TObjectDoubleHashMap<TurboID_Channel_Ori>();
						for (final TurboIDExperimentType expType : TurboIDExperimentType.values()) {
							originalIntensitiesFromTurboID.putAll(protein.getOriginalIntensities(expType));
						}
						// for array1, using row
						final TurboID_Channel_Ori channel1 = TurboID_Channel_Ori.getByName(channels.get(row));
						if (channel1 != null) {
							val1 = originalIntensitiesFromTurboID.get(channel1);
							if (applyLog2) {
								val1 = Maths.log(val1, 2);
							}

						} else {
							// if it is not returned by valueOf is because is
							// the TMT8_Nu
							final String acc = protein.getAcc();
							final ProteinFromIP proteinFromIP = ipExperiment.get(acc);
							if (proteinFromIP != null) {
								val1 = proteinFromIP.getNormalizedIntensity(IPExperimentType.NU);
								if (applyLog2) {
									val1 = Maths.log(val1, 2);
								}
							}

						}

						// for array2, using col
						final TurboID_Channel_Ori channel2 = TurboID_Channel_Ori.getByName(channels.get(col));
						if (channel2 != null) {
							val2 = originalIntensitiesFromTurboID.get(channel2);
							if (applyLog2) {
								val2 = Maths.log(val2, 2);
							}

						} else {
							// if it is not returned by valueOf is because is
							// the TMT8_Nu
							final String acc = protein.getAcc();
							final ProteinFromIP proteinFromIP = ipExperiment.get(acc);
							if (proteinFromIP != null) {
								val2 = proteinFromIP.getNormalizedIntensity(IPExperimentType.NU);
								if (applyLog2) {
									val2 = Maths.log(val2, 2);
								}
							}

						}
						// look if is significant in any bait in both
						// experiments
						if (significantlySpecificOnly) {
							if (channel1 != null) {
								final boolean specificInchannel1 = protein.isSpecific(channel1.getExpType(),
										channel1.getReplicate(),
										turboIDExperiment.getLog2PseudoSPCThreshold(channel1.getReplicate()));
								if (!specificInchannel1) {
									continue;
								}
							}
							if (channel2 != null) {
								final boolean specificInchannel2 = protein.isSpecific(channel2.getExpType(),
										channel2.getReplicate(),
										turboIDExperiment.getLog2PseudoSPCThreshold(channel2.getReplicate()));
								if (!specificInchannel2) {
									continue;
								}
							}
						}
					}
					if (val1 != null && val2 != null && !Double.isNaN(val1) && !Double.isNaN(val2)
							&& !Double.isInfinite(val1) && !Double.isInfinite(val2)) {
						array1.add(val1);
						array2.add(val2);
					}
				}
				final double correlationCoefficient = Maths.correlationCoefficient(array1, array2);
				log.info("Correlation " + channels.get(row) + " - " + channels.get(col) + " =\t"
						+ correlationCoefficient);
				matrix.set(row, col, correlationCoefficient);
			}
		}

		final FileWriter fw = new FileWriter(outputFile);

		for (final String channel : channels) {
			fw.write("\t" + channel);
		}
		fw.write("\n");
		for (row = 0; row < channels.size(); row++) {
			fw.write(channels.get(row) + "\t");
			for (col = 0; col < channels.size(); col++) {
				final double d = matrix.get(row, col);
				if (!Double.isNaN(d)) {
					fw.write(String.valueOf(d));
				}
				fw.write("\t ");
			}
			fw.write("\n");
		}
		fw.close();
		log.info("File written at: " + outputFile.getAbsolutePath());
	}

	private void exportWholeNucleiFractionFromIP_vs_TurboIDOnlyFromTurboID(TurboIDExperiment turboIDExperiment,
			IPExperiment ipExperiment, File outputFile2, boolean applyLog2, boolean useNormalizedIntensity,
			boolean onlyTM, boolean onlyNonTM, boolean significantlySpecificOnly) throws IOException {
		final FileWriter fw = new FileWriter(outputFile2);
		fw.write("accession" + "\t" + "gene" + "\t" + "transmembrane" + "\t" + "TMT8_Nu" + "\t" + "A1" + "\t"
				+ "TMT8_Nu" + "\t" + "A2" + "\t" + "TMT8_Nu" + "\t" + "A3" + "\t" + "TMT8_Nu" + "\t" + "A4" + "\t"
				+ "TMT8_Nu" + "\t" + "B1" + "\t" + "TMT8_Nu" + "\t" + "B2" + "\t" + "TMT8_Nu" + "\t" + "B3" + "\t"
				+ "TMT8_Nu" + "\t" + "B4" + "\t" + "TMT8_Nu" + "\t" + "C1" + "\t" + "TMT8_Nu" + "\t" + "C2" + "\t"
				+ "TMT8_Nu" + "\t" + "C3" + "\t" + "TMT8_Nu" + "\t" + "C4" + "\t" + "TMT8_Nu" + "\t" + "D1" + "\t"
				+ "TMT8_Nu" + "\t" + "D2" + "\t" + "TMT8_Nu" + "\t" + "D3" + "\t" + "TMT8_Nu" + "\t" + "D4" + "\t"
				+ "TMT8_Nu" + "\t" + "E1" + "\t" + "TMT8_Nu" + "\t" + "E2" + "\t" + "TMT8_Nu" + "\t" + "E3" + "\t"
				+ "TMT8_Nu" + "\t" + "E4" + "\n");
		final Collection<ProteinFromTurboID> values = turboIDExperiment.values();
		for (final ProteinFromTurboID protein : values) {
			if (onlyNonTM && protein.isTransmembrane()) {
				continue;
			}
			if (onlyTM && !protein.isTransmembrane()) {
				continue;
			}
			final String acc = protein.getAcc();
			final ProteinFromIP proteinFromIP = ipExperiment.get(acc);
			if (proteinFromIP != null) {
				Double intensity = proteinFromIP.getNormalizedIntensity(IPExperimentType.NU);
				if (!useNormalizedIntensity) {
					intensity = proteinFromIP.getOriginalIntensity(IPExperimentType.NU);
				}
				if (applyLog2) {
					intensity = Maths.log(intensity, 2);
				}
				//

				if (useNormalizedIntensity) {
					final TObjectDoubleHashMap<TurboID_Channel_Norm> intensitiesFromTurboID = new TObjectDoubleHashMap<TurboID_Channel_Norm>();

					for (final TurboIDExperimentType expType : TurboIDExperimentType.values()) {
						intensitiesFromTurboID.putAll(protein.getNormalizedIntensities(expType));
					}

					// A1-4
					String a1 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Norm.nmT2_A1), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Norm.nmT2_A1.getExpType(),
								TurboID_Channel_Norm.nmT2_A1.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Norm.nmT2_A1.getReplicate()))) {
							a1 = "";
						}

					}
					String a2 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Norm.nmT3_A2), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Norm.nmT3_A2.getExpType(),
								TurboID_Channel_Norm.nmT3_A2.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Norm.nmT3_A2.getReplicate()))) {
							a2 = "";
						}

					}
					//
					String a3 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Norm.nmT2_A3), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Norm.nmT2_A3.getExpType(),
								TurboID_Channel_Norm.nmT2_A3.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Norm.nmT2_A3.getReplicate()))) {
							a3 = "";
						}

					}
					//
					String a4 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Norm.nmT3_A4), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Norm.nmT3_A4.getExpType(),
								TurboID_Channel_Norm.nmT3_A4.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Norm.nmT3_A4.getReplicate()))) {
							a4 = "";
						}

					}
					//
					// B1-4
					String b1 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Norm.nmT4_B1), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Norm.nmT4_B1.getExpType(),
								TurboID_Channel_Norm.nmT4_B1.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Norm.nmT4_B1.getReplicate()))) {
							b1 = "";
						}

					}
					//
					String b2 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Norm.nmT5_B2), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Norm.nmT5_B2.getExpType(),
								TurboID_Channel_Norm.nmT5_B2.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Norm.nmT5_B2.getReplicate()))) {
							b2 = "";
						}

					}
					//
					String b3 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Norm.nmT4_B3), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Norm.nmT4_B3.getExpType(),
								TurboID_Channel_Norm.nmT4_B3.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Norm.nmT4_B3.getReplicate()))) {
							b3 = "";
						}

					}
					//
					String b4 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Norm.nmT5_B4), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Norm.nmT5_B4.getExpType(),
								TurboID_Channel_Norm.nmT5_B4.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Norm.nmT5_B4.getReplicate()))) {
							b4 = "";
						}

					}
					//
					// C1-4
					String c1 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Norm.nmT6_C1), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Norm.nmT6_C1.getExpType(),
								TurboID_Channel_Norm.nmT6_C1.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Norm.nmT6_C1.getReplicate()))) {
							c1 = "";
						}

					}
					//
					String c2 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Norm.nmT7_C2), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Norm.nmT7_C2.getExpType(),
								TurboID_Channel_Norm.nmT7_C2.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Norm.nmT7_C2.getReplicate()))) {
							c2 = "";
						}

					}
					//
					String c3 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Norm.nmT6_C3), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Norm.nmT6_C3.getExpType(),
								TurboID_Channel_Norm.nmT6_C3.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Norm.nmT6_C3.getReplicate()))) {
							c3 = "";
						}

					}
					//
					String c4 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Norm.nmT7_C4), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Norm.nmT7_C4.getExpType(),
								TurboID_Channel_Norm.nmT7_C4.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Norm.nmT7_C4.getReplicate()))) {
							c4 = "";
						}

					}
					//
					// D1-4
					String d1 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Norm.nmT8_D1), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Norm.nmT8_D1.getExpType(),
								TurboID_Channel_Norm.nmT8_D1.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Norm.nmT8_D1.getReplicate()))) {
							d1 = "";
						}

					}
					//
					String d2 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Norm.nmT9_D2), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Norm.nmT9_D2.getExpType(),
								TurboID_Channel_Norm.nmT9_D2.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Norm.nmT9_D2.getReplicate()))) {
							d2 = "";
						}

					}
					//
					String d3 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Norm.nmT8_D3), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Norm.nmT8_D3.getExpType(),
								TurboID_Channel_Norm.nmT8_D3.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Norm.nmT8_D3.getReplicate()))) {
							d3 = "";
						}

					}
					//
					String d4 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Norm.nmT9_D4), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Norm.nmT9_D4.getExpType(),
								TurboID_Channel_Norm.nmT9_D4.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Norm.nmT9_D4.getReplicate()))) {
							d4 = "";
						}

					}
					//
					// E1-4
					String e1 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Norm.nmT10_E1), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Norm.nmT10_E1.getExpType(),
								TurboID_Channel_Norm.nmT10_E1.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Norm.nmT10_E1.getReplicate()))) {
							e1 = "";
						}

					}
					//
					String e2 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Norm.nmT11_E2), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Norm.nmT11_E2.getExpType(),
								TurboID_Channel_Norm.nmT11_E2.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Norm.nmT11_E2.getReplicate()))) {
							e2 = "";
						}

					}
					//
					String e3 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Norm.nmT10_E3), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Norm.nmT10_E3.getExpType(),
								TurboID_Channel_Norm.nmT10_E3.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Norm.nmT10_E3.getReplicate()))) {
							e3 = "";
						}

					}
					//
					String e4 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Norm.nmT11_E4), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Norm.nmT11_E4.getExpType(),
								TurboID_Channel_Norm.nmT11_E4.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Norm.nmT11_E4.getReplicate()))) {
							e4 = "";
						}

					}
					//
					fw.write(acc + "\t" + proteinFromIP.getGene() + "\t" + proteinFromIP.isTransmembrane() + "\t"
							+ intensity + "\t" + a1 + "\t" + a2 + "\t" + a3 + "\t" + a4 + "\t" + b1 + "\t" + b2 + "\t"
							+ b3 + "\t" + b4 + "\t" + c1 + "\t" + c2 + "\t" + c3 + "\t" + c4 + "\t" + d1 + "\t" + d2
							+ "\t" + d3 + "\t" + d4 + "\t" + e1 + "\t" + e2 + "\t" + e3 + "\t" + e4 + "\n");
				} else {
					final TObjectDoubleHashMap<TurboID_Channel_Ori> intensitiesFromTurboID = new TObjectDoubleHashMap<TurboID_Channel_Ori>();

					for (final TurboIDExperimentType expType : TurboIDExperimentType.values()) {

						intensitiesFromTurboID.putAll(protein.getOriginalIntensities(expType));

					}

					// use of original intensities
					// A1-4
					String a1 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Ori.T2_A1), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Ori.T2_A1.getExpType(),
								TurboID_Channel_Ori.T2_A1.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Ori.T2_A1.getReplicate()))) {
							a1 = "";
						}

					}
					String a2 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Ori.T3_A2), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Ori.T3_A2.getExpType(),
								TurboID_Channel_Ori.T3_A2.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Ori.T3_A2.getReplicate()))) {
							a2 = "";
						}

					}
					//
					String a3 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Ori.T2_A3), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Ori.T2_A3.getExpType(),
								TurboID_Channel_Ori.T2_A3.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Ori.T2_A3.getReplicate()))) {
							a3 = "";
						}

					}
					//
					String a4 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Ori.T3_A4), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Ori.T3_A4.getExpType(),
								TurboID_Channel_Ori.T3_A4.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Ori.T3_A4.getReplicate()))) {
							a4 = "";
						}

					}
					//
					// B1-4
					String b1 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Ori.T4_B1), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Ori.T4_B1.getExpType(),
								TurboID_Channel_Ori.T4_B1.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Ori.T4_B1.getReplicate()))) {
							b1 = "";
						}

					}
					//
					String b2 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Ori.T5_B2), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Ori.T5_B2.getExpType(),
								TurboID_Channel_Ori.T5_B2.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Ori.T5_B2.getReplicate()))) {
							b2 = "";
						}

					}
					//
					String b3 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Ori.T4_B3), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Ori.T4_B3.getExpType(),
								TurboID_Channel_Ori.T4_B3.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Ori.T4_B3.getReplicate()))) {
							b3 = "";
						}

					}
					//
					String b4 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Ori.T5_B4), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Ori.T5_B4.getExpType(),
								TurboID_Channel_Ori.T5_B4.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Ori.T5_B4.getReplicate()))) {
							b4 = "";
						}

					}
					//
					// C1-4
					String c1 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Ori.T6_C1), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Ori.T6_C1.getExpType(),
								TurboID_Channel_Ori.T6_C1.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Ori.T6_C1.getReplicate()))) {
							c1 = "";
						}

					}
					//
					String c2 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Ori.T7_C2), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Ori.T7_C2.getExpType(),
								TurboID_Channel_Ori.T7_C2.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Ori.T7_C2.getReplicate()))) {
							c2 = "";
						}

					}
					//
					String c3 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Ori.T6_C3), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Ori.T6_C3.getExpType(),
								TurboID_Channel_Ori.T6_C3.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Ori.T6_C3.getReplicate()))) {
							c3 = "";
						}

					}
					//
					String c4 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Ori.T7_C4), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Ori.T7_C4.getExpType(),
								TurboID_Channel_Ori.T7_C4.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Ori.T7_C4.getReplicate()))) {
							c4 = "";
						}

					}
					//
					// D1-4
					String d1 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Ori.T8_D1), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Ori.T8_D1.getExpType(),
								TurboID_Channel_Ori.T8_D1.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Ori.T8_D1.getReplicate()))) {
							d1 = "";
						}

					}
					//
					String d2 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Ori.T9_D2), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Ori.T9_D2.getExpType(),
								TurboID_Channel_Ori.T9_D2.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Ori.T9_D2.getReplicate()))) {
							d2 = "";
						}

					}
					//
					String d3 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Ori.T8_D3), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Ori.T8_D3.getExpType(),
								TurboID_Channel_Ori.T8_D3.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Ori.T8_D3.getReplicate()))) {
							d3 = "";
						}

					}
					//
					String d4 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Ori.T9_D4), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Ori.T9_D4.getExpType(),
								TurboID_Channel_Ori.T9_D4.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Ori.T9_D4.getReplicate()))) {
							d4 = "";
						}

					}
					//
					// E1-4
					String e1 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Ori.T10_E1), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Ori.T10_E1.getExpType(),
								TurboID_Channel_Ori.T10_E1.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Ori.T10_E1.getReplicate()))) {
							e1 = "";
						}

					}
					//
					String e2 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Ori.T11_E2), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Ori.T11_E2.getExpType(),
								TurboID_Channel_Ori.T11_E2.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Ori.T11_E2.getReplicate()))) {
							e2 = "";
						}

					}
					//
					String e3 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Ori.T10_E3), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Ori.T10_E3.getExpType(),
								TurboID_Channel_Ori.T10_E3.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Ori.T10_E3.getReplicate()))) {
							e3 = "";
						}

					}
					//
					String e4 = parseDouble(intensitiesFromTurboID.get(TurboID_Channel_Ori.T11_E4), intensity,
							applyLog2);
					if (significantlySpecificOnly) {
						if (!protein.isSpecific(TurboID_Channel_Ori.T11_E4.getExpType(),
								TurboID_Channel_Ori.T11_E4.getReplicate(), turboIDExperiment
										.getLog2PseudoSPCThreshold(TurboID_Channel_Ori.T11_E4.getReplicate()))) {
							e4 = "";
						}

					}
					//
					fw.write(acc + "\t" + proteinFromIP.getGene() + "\t" + proteinFromIP.isTransmembrane() + "\t"
							+ intensity + "\t" + a1 + "\t" + a2 + "\t" + a3 + "\t" + a4 + "\t" + b1 + "\t" + b2 + "\t"
							+ b3 + "\t" + b4 + "\t" + c1 + "\t" + c2 + "\t" + c3 + "\t" + c4 + "\t" + d1 + "\t" + d2
							+ "\t" + d3 + "\t" + d4 + "\t" + e1 + "\t" + e2 + "\t" + e3 + "\t" + e4 + "\n");
				}
			}
		}
		fw.close();
		log.info("File written at: " + outputFile2.getAbsolutePath());
	}

	private String parseDouble(double d, Double normalizedIntensity, boolean applyLog2) {
		if (Double.isNaN(d) || Double.compare(d, 0.0) == 0) {
			return "\t";
		}
		if (applyLog2) {
			return Maths.log(d, 2) + "\t" + normalizedIntensity;
		}
		return d + "\t" + normalizedIntensity;
	}

	private IPExperiment createIPExperiment() throws IOException {
		final Map<String, String> genesByACC = new THashMap<String, String>();
		final IPExperiment ipExperiment = new IPExperiment();
		final ExcelReader reader = new ExcelReader(excelFile_TMT8_EMD_IP, 2, 3);
		final Set<String> genes = new THashSet<String>();
		final int sheetIndex = reader.getWorkbook().getSheetIndex(sheet_TMT8_EMD_IP);
		final int columnIndexForACC = reader.getColumnIndex(sheetIndex, ACCESSION);
		final int columnIndexForDescription = reader.getColumnIndex(sheetIndex, DESCRIPTION);
		final int columnIndexForSPC = reader.getColumnIndex(sheetIndex, PSM);
		int numRow = 1;
		while (true) {
			final String acc = reader.getStringValue(sheetIndex, numRow, columnIndexForACC);
			if (acc == null) {
				break;
			}
			final String description = reader.getStringValue(sheetIndex, numRow, columnIndexForDescription);
			final int spc = Double.valueOf(reader.getStringValue(sheetIndex, numRow, columnIndexForSPC)).intValue();
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
			final String transmembraneString = reader.getStringValue(sheetIndex, numRow, columnIndex);
			if (transmembraneString != null && !"".equals(transmembraneString)) {
				isTransmembrane = Boolean.valueOf(transmembraneString);
			}
			ProteinFromIP protein = null;
			if (ipExperiment.containsKey(acc)) {
				protein = ipExperiment.get(acc);
			} else {
				protein = new ProteinFromIP(acc, gene, isTransmembrane);
				ipExperiment.put(protein.getAcc(), protein);
			}
			protein.setSpc(spc);
			numRow++;
		}

		log.info(ipExperiment.size() + " proteins in total");

		numRow = 1;
		while (true) {
			final String acc = reader.getStringValue(sheetIndex, numRow, columnIndexForACC);
			if (acc == null) {
				break;
			}

			final ProteinFromIP protein = ipExperiment.get(acc);

			final Set<IP_Channel_Norm> notFound = new THashSet<IP_Channel_Norm>();
			int counter = 0;
			for (final IP_Channel_Norm channel : IP_Channel_Norm.values()) {
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
				protein.addNormalizedIntensity(intensity, channel);
			}
			if (counter != IP_Channel_Norm.values().length) {
				for (final IP_Channel_Norm string : notFound) {
					log.info(string);
				}
				throw new IllegalArgumentException("Some columns were not found");
			}
			///////////////////////////////
			final Set<IP_Channel_Ori> notFound2 = new THashSet<IP_Channel_Ori>();
			counter = 0;
			for (final IP_Channel_Ori channel : IP_Channel_Ori.values()) {
				final int columnIndex = reader.getColumnIndex(sheetIndex, channel.name());
				if (columnIndex == -1) {
					notFound2.add(channel);
					continue;
				}
				counter++;
				final String stringValue = reader.getNumberValue(sheetIndex, numRow, columnIndex);
				Double intensity = Double.NaN;
				if (stringValue != null && !stringValue.equals("NA")) {
					intensity = Double.valueOf(stringValue);
				}
				protein.addOriginalIntensity(intensity, channel);
			}
			if (counter != IP_Channel_Ori.values().length) {
				for (final IP_Channel_Ori string : notFound2) {
					log.info(string);
				}
				throw new IllegalArgumentException("Some columns were not found");
			}
			numRow++;
		}
		return ipExperiment;
	}

	private TurboIDExperiment readTurboIDExperiment(File file, String[] sheets) throws IOException {
		final Map<String, String> genesByACC = new THashMap<String, String>();
		final TurboIDExperiment turboIDExperiment = new TurboIDExperiment(fraction);
		final ExcelReader reader = new ExcelReader(file, 1, 4);
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
			int numRow = 2;
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
							intensity = Double.valueOf(stringValue);
						}
						final TurboIDExperimentType turboIDExperimenType = replicateID.getExpType();
						protein.addOriginalIntensity(intensity, replicateID, turboIDExperimenType);
					}
					// if (counter != TurboID_Channel_Ori.values().length / 3) {
					// for (final TurboID_Channel_Ori string : notFound2) {
					// log.info(string);
					// }
					// throw new IllegalArgumentException("Some columns were not
					// found");
					// }
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

	private void loadUniprotAnnotations(Collection<String> accs) {
		log.info("Loading Uniprot annotations for " + accs.size() + " proteins");
		final UniprotProteinRetriever upr = new UniprotProteinRetriever(null, uniprotReleasesFolder, true);
		upr.getAnnotatedProteins(accs);
		log.info("annotations loaded");
	}

	private String getOutputFileName(String prefix, boolean onlyTM, boolean onlyNonTM, boolean ignoreMissingValues,
			boolean useOfNormalizedIntensities, boolean normalizeToAverage, boolean normalizeToTurboIDChannelAverage,
			boolean log2Applied) {
		final String v5 = "";
		// if (normalizingByV5Tag) {
		// v5 = "_V5";
		// }
		final String mixChannel = "";
		// if (normalizingByMixChannel) {
		// mixChannel = "_mix";
		// }
		String tmString = "_all";
		if (onlyTM || onlyNonTM) {
			if (onlyTM) {
				tmString = "_TM";
			}
			if (onlyNonTM) {
				tmString = "_nonTM";
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
		replicateBySheets.put(sheetsTMT11Nu[1], Replicate.B);
		replicateBySheets.put(sheetsTMT11Nu[2], Replicate.Brerun);

		replicateBySheets.put(sheetsTMT11Cy[0], Replicate.CyA);
		replicateBySheets.put(sheetsTMT11Cy[1], Replicate.CyB);
	}

	private String parseForExcelConflictingGenes(String gene) {
		if (gene == null) {
			gene = "N/A";
		}
		if (gene.equals("March5")) {
			gene = "March5_";
		}
		return gene;
	}
}
