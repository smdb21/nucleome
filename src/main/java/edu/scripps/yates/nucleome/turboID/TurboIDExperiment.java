package edu.scripps.yates.nucleome.turboID;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.collections.comparators.ReverseComparator;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.apache.log4j.Logger;

import edu.scripps.yates.utilities.maths.Maths;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;

public class TurboIDExperiment extends THashMap<String, ProteinFromTurboID> {
	private final static Logger log = Logger.getLogger(TurboIDExperiment.class);

	private static final String V5TAG = "V5tag";
	private ProteinFromTurboID v5TagProtein;

	private final TObjectDoubleHashMap<TurboID_Channel_Ori> sumIntensitiesByChannel = new TObjectDoubleHashMap<TurboID_Channel_Ori>();

	private final TObjectDoubleHashMap<Replicate> log2PseudoSPCThresholdPerReplicate = new TObjectDoubleHashMap<Replicate>();
	private final THashMap<Replicate, Gaussian> gaussianFitPerReplicate = new THashMap<Replicate, Gaussian>();
	private final THashMap<Replicate, List<WeightedObservedPoint>> experimentalHistogramPerReplicate = new THashMap<Replicate, List<WeightedObservedPoint>>();

	private Gaussian gaussianFit;

	private List<WeightedObservedPoint> experimentalHistogram;

	private double log2PseudoSPCThreshold;

	private final TurboIDFraction fraction;

	// private boolean normalizingByMixChannel = false;
	// private boolean normalizingByV5Tag = false;;

	public TurboIDExperiment(TurboIDFraction fraction) {
		this.fraction = fraction;
	}

	public void normalizeByV5Tag(boolean applyLog2) {
		if (v5TagProtein == null) {
			throw new IllegalArgumentException("V5Tag was not found in input data");
		}
		log.info("Normalizing by V5 tag");
		// normalizingByV5Tag = true;
		normalizeNormalizedByV5Tag(applyLog2);
		normalizeOriginalByV5Tag(applyLog2);
	}

	public TurboIDFraction getFraction() {
		return this.fraction;
	}

	private void normalizeOriginalByV5Tag(boolean applyLog2) {
		System.out.println("Normalizing original intensities by the V5 tag");
		for (final TurboIDExperimentType type : TurboIDExperimentType.values()) {
			final TObjectDoubleHashMap<TurboID_Channel_Ori> v5Intensities = new TObjectDoubleHashMap<TurboID_Channel_Ori>();
			v5Intensities.putAll(v5TagProtein.getOriginalIntensities(type));
			for (final ProteinFromTurboID protein : values()) {
				final TObjectDoubleHashMap<TurboID_Channel_Ori> nonNormalizedIntensities = protein
						.getOriginalIntensities(type);
				for (final TurboID_Channel_Ori replicateID : v5Intensities.keySet()) {
					if (nonNormalizedIntensities.containsKey(replicateID)) {
						final double intensity = nonNormalizedIntensities.get(replicateID);
						if (!Double.isNaN(intensity)) {
							double normalizedIntensity = intensity / v5Intensities.get(replicateID);
							if (applyLog2) {
								normalizedIntensity = Maths.log(normalizedIntensity, 2);
							}
							nonNormalizedIntensities.put(replicateID, normalizedIntensity);
						}
					}
				}
			}
			// print them
			for (final TurboID_Channel_Ori channel : TurboID_Channel_Ori.values()) {
				if (v5Intensities.containsKey(channel)) {
					System.out.println(
							"V5 tag: " + type.name() + "\t" + channel.name() + "\t" + v5Intensities.get(channel));
				}
			}
		}
	}

	private void normalizeNormalizedByV5Tag(boolean applyLog2) {
		System.out.println("Normalizing IP2-normalized intensities by the V5 tag");
		for (final TurboIDExperimentType type : TurboIDExperimentType.values()) {
			final TObjectDoubleHashMap<TurboID_Channel_Norm> v5Intensities = new TObjectDoubleHashMap<TurboID_Channel_Norm>();
			v5Intensities.putAll(v5TagProtein.getNormalizedIntensities(type));
			for (final ProteinFromTurboID protein : values()) {
				final TObjectDoubleHashMap<TurboID_Channel_Norm> nonNormalizedIntensities = protein
						.getNormalizedIntensities(type);
				for (final TurboID_Channel_Norm replicateID : v5Intensities.keySet()) {
					if (nonNormalizedIntensities.containsKey(replicateID)) {
						final double intensity = nonNormalizedIntensities.get(replicateID);
						if (!Double.isNaN(intensity)) {
							double normalizedIntensity = intensity / v5Intensities.get(replicateID);
							if (applyLog2) {
								normalizedIntensity = Maths.log(normalizedIntensity, 2);
							}
							nonNormalizedIntensities.put(replicateID, normalizedIntensity);
						}
					}
				}
			}
			// print them
			for (final TurboID_Channel_Norm channel : TurboID_Channel_Norm.values()) {
				if (v5Intensities.containsKey(channel)) {
					System.out.println(
							"V5 tag: " + type.name() + "\t" + channel.name() + "\t" + v5Intensities.get(channel));
				}
			}
		}

	}

	@Override
	public ProteinFromTurboID put(String acc, ProteinFromTurboID protein) {
		if (acc.toLowerCase().contains(V5TAG.toLowerCase())) {
			v5TagProtein = protein;
			v5TagProtein.setV5Tag(true);
		}
		return super.put(acc, protein);
	}

	public void normalizeByMixChannel() {
		log.info("Normalizing by mix channel");
		// normalizingByMixChannel = true;
		normalizeNormalizedByMixChannel();
		normalizeOriginalByMixChannel();

	}

	private void normalizeOriginalByMixChannel() {
		for (final ProteinFromTurboID protein : values()) {
			final TObjectDoubleHashMap<TurboID_Channel_Ori> mixIntensities = protein
					.getOriginalIntensities(TurboIDExperimentType.MIX);

			for (final TurboIDExperimentType type : TurboIDExperimentType.values()) {
				if (type == TurboIDExperimentType.MIX) {
					continue;
				}
				final TObjectDoubleHashMap<TurboID_Channel_Ori> nonNormalizedIntensities = protein
						.getOriginalIntensities(type);
				for (final TurboID_Channel_Ori replicateID : nonNormalizedIntensities.keySet()) {
					final Replicate replicate = replicateID.getReplicate();

					boolean found = false;
					for (final TurboID_Channel_Ori mixIntensityChannel : mixIntensities.keySet()) {
						if (mixIntensityChannel.getReplicate() == replicate) {
							found = true;
							final double mixIntensity = mixIntensities.get(mixIntensityChannel);
							final double intensity = nonNormalizedIntensities.get(replicateID);
							if (!Double.isNaN(intensity)) {
								final double normalizedIntensity = intensity / mixIntensity;
								nonNormalizedIntensities.put(replicateID, normalizedIntensity);
							}
						}
					}
					if (!found) {
						// set all to N/A because that is because there is no
						// mix channel
						nonNormalizedIntensities.put(replicateID, Double.NaN);
					}
				}
			}
		}
	}

	private void normalizeNormalizedByMixChannel() {
		for (final ProteinFromTurboID protein : values()) {
			final TObjectDoubleHashMap<TurboID_Channel_Norm> mixIntensities = protein
					.getNormalizedIntensities(TurboIDExperimentType.MIX);

			for (final TurboIDExperimentType type : TurboIDExperimentType.values()) {
				if (type == TurboIDExperimentType.MIX) {
					continue;
				}
				final TObjectDoubleHashMap<TurboID_Channel_Norm> nonNormalizedIntensities = protein
						.getNormalizedIntensities(type);
				for (final TurboID_Channel_Norm replicateID : nonNormalizedIntensities.keySet()) {
					final Replicate replicate = replicateID.getReplicate();

					boolean found = false;
					for (final TurboID_Channel_Norm mixIntensityChannel : mixIntensities.keySet()) {
						if (mixIntensityChannel.getReplicate() == replicate) {
							found = true;
							final double mixIntensity = mixIntensities.get(mixIntensityChannel);
							final double intensity = nonNormalizedIntensities.get(replicateID);
							if (!Double.isNaN(intensity)) {
								final double normalizedIntensity = intensity / mixIntensity;
								nonNormalizedIntensities.put(replicateID, normalizedIntensity);
							}
						}
					}
					if (!found) {
						// set all to N/A because that is because there is no
						// mix channel
						nonNormalizedIntensities.put(replicateID, Double.NaN);
					}
				}
			}
		}
	}

	public void normalizeNormalizedByAvgOfBaitChannels(boolean ignoreMissingValues, boolean applyLog2) {
		for (final ProteinFromTurboID protein : values()) {
			final TDoubleArrayList toAverage = new TDoubleArrayList();
			for (final TurboIDExperimentType turboIDExperimentType : TurboIDExperimentType.values()) {
				if (!turboIDExperimentType.isBait()) {
					continue; // not use mix channel
				}
				final double intensity = protein.getAverageFromNormalized(turboIDExperimentType, ignoreMissingValues);
				if (!Double.isNaN(intensity)) {
					toAverage.add(intensity);
				}
			}
			final double averageOfBaitChannels = Maths.mean(toAverage);

			for (final TurboIDExperimentType type : TurboIDExperimentType.values()) {
				if (!type.isBait()) {
					continue;
				}
				final TObjectDoubleHashMap<TurboID_Channel_Norm> nonNormalizedIntensities = protein
						.getNormalizedIntensities(type);
				for (final TurboID_Channel_Norm replicateID : nonNormalizedIntensities.keySet()) {
					double value = nonNormalizedIntensities.get(replicateID) / averageOfBaitChannels;
					if (applyLog2) {
						value = Maths.log(value, 2);
					}
					nonNormalizedIntensities.put(replicateID, value);
				}
			}
		}
	}

	public void normalizeToTurboIDChannelAverage(boolean ignoreMissingValues, boolean applyLog2,
			Collection<Replicate> replicatesToUse) {
		for (final ProteinFromTurboID protein : values()) {

			final double turboIDOnly = protein.getAverageFromNormalized(TurboIDExperimentType.TURBOID_ONLY,
					ignoreMissingValues, replicatesToUse);

			for (final TurboIDExperimentType type : TurboIDExperimentType.values()) {
				if (!type.isBait()) {
					continue;
				}
				final TObjectDoubleHashMap<TurboID_Channel_Norm> nonNormalizedIntensities = protein
						.getNormalizedIntensities(type);
				for (final TurboID_Channel_Norm replicateID : nonNormalizedIntensities.keySet()) {
					double value = nonNormalizedIntensities.get(replicateID) / turboIDOnly;
					if (applyLog2) {
						value = Maths.log(value, 2);
					}
					nonNormalizedIntensities.put(replicateID, value);
				}
			}
		}
	}

	// public boolean isNormalizingByMixChannel() {
	// return normalizingByMixChannel;
	// }

	// public boolean isNormalizingByV5Tag() {
	// return normalizingByV5Tag;
	// }

	public File exportToFileForClustering(File outputFile, boolean onlyTM, boolean onlyNonTM,
			boolean ignoreMissingValues, boolean useNormalizedIntensities) {

		// applyLogs();
		dealWithInfinities();

		FileWriter fw = null;
		try {
			fw = new FileWriter(outputFile);
			// header
			fw.write("ACC\t");
			fw.write("Gene\t");
			fw.write("Transmembrane\t");
			// spc
			for (final Replicate replicate : Replicate.values(fraction)) {
				fw.write("SPC " + replicate.name() + "\t");
			}
			for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaits()) {
				fw.write("specific_" + bait.name() + "\t");
			}
			StringBuilder sb = new StringBuilder();
			for (final TurboIDExperimentType turboIDExperimentType : TurboIDExperimentType.values()) {
				if (turboIDExperimentType != TurboIDExperimentType.MIX) {
					if (!"".equals(sb.toString())) {
						sb.append("\t");
					}
					sb.append(turboIDExperimentType);
					if (turboIDExperimentType == TurboIDExperimentType.TURBOID_ONLY) {

						for (final Replicate replicate : Replicate.values(fraction)) {
							for (final TurboID_Channel_Norm channel : TurboID_Channel_Norm.values()) {
								if (channel.getReplicate() == replicate) {
									if (channel.name().contains("10") || channel.name().contains("11")) {

										if (!"".equals(sb.toString())) {
											sb.append("\t");
										}
										sb.append(channel.getReplicate() + "-" + channel.name());
									}
								}
							}

						}
					}
				}
			}
			fw.write(sb.toString());
			fw.write("\n");
			for (final ProteinFromTurboID protein : values()) {
				if (onlyTM && !protein.isTransmembrane()) {
					continue;
				}
				if (onlyNonTM && protein.isTransmembrane()) {
					continue;
				}
				try {
					fw.write(protein.getAcc() + "\t");
					fw.write(protein.getGene() + "\t");
					fw.write(protein.isTransmembrane() + "\t");
					// spc
					for (final Replicate replicate : Replicate.values(fraction)) {
						fw.write(protein.getSpc(replicate) + "\t");
					}
					for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaits()) {
						fw.write(protein.isSpecificInBothExperiments(bait, getLog2PseudoSPCThreshold(Replicate.A),
								getLog2PseudoSPCThreshold(Replicate.B)) + "\t");
					}
					sb = new StringBuilder();
					for (final TurboIDExperimentType turboIDExperimentType : TurboIDExperimentType.values()) {
						if (turboIDExperimentType != TurboIDExperimentType.MIX) {
							double avg = Double.NaN;
							if (useNormalizedIntensities) {
								avg = protein.getAverageFromNormalized(turboIDExperimentType, ignoreMissingValues);
							} else {
								avg = protein.getAverageFromOriginal(turboIDExperimentType, ignoreMissingValues);
							}
							if (!"".equals(sb.toString())) {
								sb.append("\t");
							}
							sb.append(avg);
							if (turboIDExperimentType == TurboIDExperimentType.TURBOID_ONLY) {
								final TObjectDoubleHashMap<TurboID_Channel_Norm> normalizedIntensities = protein
										.getNormalizedIntensities(turboIDExperimentType);
								for (final Replicate replicate : Replicate.values(fraction)) {
									for (final TurboID_Channel_Norm channel : TurboID_Channel_Norm.values()) {
										if (channel.getReplicate() == replicate) {
											if (channel.name().contains("10") || channel.name().contains("11")) {
												final double intensity = normalizedIntensities.get(channel);
												if (!"".equals(sb.toString())) {
													sb.append("\t");
												}
												sb.append(intensity);
											}
										}
									}

								}
							}
						}
					}
					fw.write(sb.toString());
				} finally {
					fw.write("\n");
				}
			}
		} catch (final IOException e) {
			e.printStackTrace();
		} finally {
			if (fw != null) {
				try {
					fw.close();
				} catch (final IOException e) {
				}
			}
		}
		return outputFile;
	}

	// private void applyLogs() {
	// for (final ProteinFromTurboID protein : values()) {
	// for (final TurboIDExperimentType expType :
	// TurboIDExperimentType.values()) {
	// final TObjectDoubleHashMap<Channel_Norm> intensities =
	// protein.getNormalizedIntensities(expType);
	// for (final Channel_Norm channel : intensities.keySet()) {
	// final double d = intensities.get(channel);
	// intensities.put(channel, Maths.log(d, 2));
	// }
	// final TObjectDoubleHashMap<Channel_Ori> intensities2 =
	// protein.getOriginalIntensities(expType);
	// for (final Channel_Ori channel : intensities2.keySet()) {
	// final double d = intensities2.get(channel);
	// intensities2.put(channel, Maths.log(d, 2));
	// }
	// }
	// }
	// }

	public void dealWithInfinities() {
		for (final TurboIDExperimentType expType : TurboIDExperimentType.values()) {
			Double maxNorm = -Double.MAX_VALUE;
			Double minNorm = Double.MAX_VALUE;
			Double maxOri = -Double.MAX_VALUE;
			Double minOri = Double.MAX_VALUE;
			for (final ProteinFromTurboID protein : values()) {
				final TObjectDoubleHashMap<TurboID_Channel_Norm> intensities = protein
						.getNormalizedIntensities(expType);
				for (final TurboID_Channel_Norm channel : intensities.keySet()) {
					final double d = intensities.get(channel);
					if (Double.isFinite(d)) {
						if (d > maxNorm) {
							maxNorm = d;
						}
						if (d < minNorm) {
							minNorm = d;
						}
					}
				}
				final TObjectDoubleHashMap<TurboID_Channel_Ori> intensities2 = protein.getOriginalIntensities(expType);
				for (final TurboID_Channel_Ori channel : intensities2.keySet()) {
					final double d = intensities2.get(channel);
					if (Double.isFinite(d)) {
						if (d > maxOri) {
							maxOri = d;
						}
						if (d < minOri) {
							minOri = d;
						}
					}
				}
			}
			log.info("Maximum from normalized for " + expType + ": " + maxNorm);
			log.info("Minimum from normalized for " + expType + ": " + minNorm);
			log.info("Maximum from original for " + expType + ": " + maxOri);
			log.info("Minimum from original for " + expType + ": " + minOri);
			//
			for (final ProteinFromTurboID protein : values()) {
				final TObjectDoubleHashMap<TurboID_Channel_Norm> intensities = protein
						.getNormalizedIntensities(expType);
				for (final TurboID_Channel_Norm channel : intensities.keySet()) {
					final double d = intensities.get(channel);
					if (Double.isInfinite(d)) {
						if (Double.compare(Double.NEGATIVE_INFINITY, d) == 0) {
							intensities.put(channel, minNorm);
						}
						if (Double.compare(Double.POSITIVE_INFINITY, d) == 0) {
							intensities.put(channel, maxNorm);
						}
					}
				}
				final TObjectDoubleHashMap<TurboID_Channel_Ori> intensities2 = protein.getOriginalIntensities(expType);
				for (final TurboID_Channel_Ori channel : intensities2.keySet()) {
					final double d = intensities2.get(channel);
					if (Double.isInfinite(d)) {
						if (Double.compare(Double.NEGATIVE_INFINITY, d) == 0) {
							intensities2.put(channel, minOri);
						}
						if (Double.compare(Double.POSITIVE_INFINITY, d) == 0) {
							intensities2.put(channel, maxOri);
						}
					}
				}
			}
		}
	}

	public static Comparator<ProteinFromTurboID> getComparatorByTotalSPC(boolean ascending) {
		final Comparator<ProteinFromTurboID> ret = new Comparator<ProteinFromTurboID>() {

			@Override
			public int compare(ProteinFromTurboID o1, ProteinFromTurboID o2) {
				return Integer.compare(o1.getSumSPCAcrossReplicates(), o2.getSumSPCAcrossReplicates());
			}
		};
		if (!ascending) {
			return new ReverseComparator(ret);
		}
		return ret;
	}

	public static Comparator<ProteinFromTurboID> getComparatorByAnova(boolean ascending, boolean useLogs,
			boolean distribSPC, boolean distribIntensity, TurboIDFraction fraction) {
		final Comparator<ProteinFromTurboID> ret = new Comparator<ProteinFromTurboID>() {

			@Override
			public int compare(ProteinFromTurboID o1, ProteinFromTurboID o2) {
				return Double.compare(o1.getAnovaPValueOverBaits(useLogs, distribSPC, distribIntensity, fraction),
						o2.getAnovaPValueOverBaits(useLogs, distribSPC, distribIntensity, fraction));
			}
		};
		if (!ascending) {
			return new ReverseComparator(ret);
		}
		return ret;
	}

	public static Comparator<ProteinFromTurboID> getComparatorBySumIntensitiesAcrossReplicates(boolean ascending,
			TurboIDExperimentType expType, boolean useNormalizedIntensities) {
		final Comparator<ProteinFromTurboID> ret = new Comparator<ProteinFromTurboID>() {

			@Override
			public int compare(ProteinFromTurboID o1, ProteinFromTurboID o2) {
				return Double.compare(o1.getSumIntesitiesAcrossReplicates(expType, useNormalizedIntensities),
						o2.getSumIntesitiesAcrossReplicates(expType, useNormalizedIntensities));
			}
		};
		if (!ascending) {
			return new ReverseComparator(ret);
		}
		return ret;
	}

	/**
	 * Gets the proteins sorted using a comparator. <br>
	 * Look for static methods that create comparators for {@link TurboIDExperiment}
	 * 
	 * @param comparator
	 * @return
	 */
	public List<ProteinFromTurboID> getProteinsSorted(Comparator<ProteinFromTurboID> comparator) {
		final List<ProteinFromTurboID> ret = new ArrayList<ProteinFromTurboID>();
		ret.addAll(values());
		Collections.sort(ret, comparator);
		return ret;
	}

	public double getSumIntensities(TurboID_Channel_Ori channelOri) {
		if (!sumIntensitiesByChannel.containsKey(channelOri)) {
			double ret = 0.0;
			for (final ProteinFromTurboID protein : values()) {
				final double intensity = protein.getOriginalIntensities(channelOri.getExpType()).get(channelOri);
				if (!Double.isNaN(intensity)) {
					ret += intensity;
				}
			}
			sumIntensitiesByChannel.put(channelOri, ret);
			return ret;
		}
		return sumIntensitiesByChannel.get(channelOri);
	}

	public void setLog2PseudoSPCThreshold(double threshold, Replicate replicate) {
		log2PseudoSPCThresholdPerReplicate.put(replicate, threshold);
	}

	public void setLog2PseudoSPCThreshold(double threshold) {
		log2PseudoSPCThreshold = threshold;
	}

	public double getLog2PseudoSPCThreshold() {
		return log2PseudoSPCThreshold;
	}

	public double getLog2PseudoSPCThreshold(Replicate replicate) {
		return log2PseudoSPCThresholdPerReplicate.get(replicate);
	}

	public void setGaussianFit(Gaussian gaussian, Replicate replicate) {
		gaussianFitPerReplicate.put(replicate, gaussian);
	}

	public Gaussian getGaussianFit(Replicate replicate) {
		return gaussianFitPerReplicate.get(replicate);
	}

	public Gaussian getGaussianFit() {
		return gaussianFit;
	}

	public void setGaussianFit(Gaussian gaussianFit) {
		this.gaussianFit = gaussianFit;
	}

	public void setExperimentalHistogram(List<WeightedObservedPoint> points, Replicate replicate) {
		experimentalHistogramPerReplicate.put(replicate, points);
	}

	public List<WeightedObservedPoint> getExperimentalHistogram(Replicate replicate) {
		return experimentalHistogramPerReplicate.get(replicate);
	}

	public List<WeightedObservedPoint> getExperimentalHistogram() {
		return experimentalHistogram;
	}

	public void setExperimentalHistogram(List<WeightedObservedPoint> experimentalHistogram) {
		this.experimentalHistogram = experimentalHistogram;
	}

}
