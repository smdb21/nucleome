package edu.scripps.yates.nucleome.turboID;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.apache.commons.math3.stat.inference.OneWayAnova;
import org.apache.commons.math3.stat.inference.TTest;
import org.apache.log4j.Logger;

import edu.scripps.yates.utilities.maths.Maths;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import gnu.trove.map.hash.TObjectIntHashMap;

public class ProteinFromTurboID {
	private final THashMap<TurboIDExperimentType, TObjectDoubleHashMap<TurboID_Channel_Norm>> normalizedIntensities = new THashMap<TurboIDExperimentType, TObjectDoubleHashMap<TurboID_Channel_Norm>>();
	private final THashMap<TurboIDExperimentType, TObjectDoubleHashMap<TurboID_Channel_Ori>> originalIntensities = new THashMap<TurboIDExperimentType, TObjectDoubleHashMap<TurboID_Channel_Ori>>();
	private final String acc;
	private final boolean isTransmembrane;
	private boolean isV5Tag = false;
	private final TObjectIntHashMap<Replicate> spcByReplicate = new TObjectIntHashMap<Replicate>();
	private final String gene;
	private int sumSPCAcrossReplicates = -1;
	private final THashMap<Replicate, TObjectDoubleHashMap<TurboID_Channel_Norm>> pseudoSpecCountsWithNormalizedIntensitiesByReplicate = new THashMap<Replicate, TObjectDoubleHashMap<TurboID_Channel_Norm>>();
	private final TObjectDoubleHashMap<TurboID_Channel_Norm> distributedIntensitiesIntensities = new TObjectDoubleHashMap<TurboID_Channel_Norm>();
	private final THashMap<Replicate, TObjectDoubleHashMap<TurboID_Channel_Ori>> pseudoSpecCountsWithOriginalIntensitiesByReplicate = new THashMap<Replicate, TObjectDoubleHashMap<TurboID_Channel_Ori>>();
	private final static Logger log = Logger.getLogger(ProteinFromTurboID.class);

	public ProteinFromTurboID(String acc, String gene, boolean isTransmembrane) {
		this.acc = acc;
		this.isTransmembrane = isTransmembrane;
		this.gene = gene;
	}

	public boolean addNormalizedIntensity(double intensity, TurboID_Channel_Norm replicateID,
			TurboIDExperimentType turboIDExperiment) {
		if (Double.isNaN(intensity)) {
			// return false;
		}
		if (Double.compare(intensity, 0.0) == 0) {
			// return false;
			intensity = Double.NaN;
		}
		if (normalizedIntensities.containsKey(turboIDExperiment)) {
			normalizedIntensities.get(turboIDExperiment).put(replicateID, intensity);
		} else {
			final TObjectDoubleHashMap<TurboID_Channel_Norm> intensitiesForThisProtein = new TObjectDoubleHashMap<TurboID_Channel_Norm>();
			intensitiesForThisProtein.put(replicateID, intensity);
			normalizedIntensities.put(turboIDExperiment, intensitiesForThisProtein);
		}
		return true;
	}

	public String getAcc() {
		return acc;
	}

	public boolean isTransmembrane() {
		return isTransmembrane;
	}

	public TObjectDoubleHashMap<TurboID_Channel_Norm> getNormalizedIntensities(TurboIDExperimentType experiment) {
		if (!normalizedIntensities.containsKey(experiment)) {
			normalizedIntensities.put(experiment, new TObjectDoubleHashMap<TurboID_Channel_Norm>());
		}
		return normalizedIntensities.get(experiment);
	}

	public TObjectDoubleHashMap<TurboID_Channel_Ori> getOriginalIntensities(TurboIDExperimentType experiment) {
		if (!originalIntensities.containsKey(experiment)) {
			originalIntensities.put(experiment, new TObjectDoubleHashMap<TurboID_Channel_Ori>());
		}
		return originalIntensities.get(experiment);
	}

	public boolean isV5Tag() {
		return isV5Tag;
	}

	public void setV5Tag(boolean isV5Tag) {
		this.isV5Tag = isV5Tag;
	}

	public double getAverageFromNormalized(TurboIDExperimentType turboIDExperimentType, boolean ignoreMissingValues,
			Collection<Replicate> replicatesToUse) {
		final TObjectDoubleHashMap<TurboID_Channel_Norm> intensities = getNormalizedIntensities(turboIDExperimentType);
		final TObjectDoubleHashMap<TurboID_Channel_Norm> toAverage = new TObjectDoubleHashMap<TurboID_Channel_Norm>();
		for (final TurboID_Channel_Norm channel : intensities.keySet()) {
			if (replicatesToUse.contains(channel.getReplicate())) {
				toAverage.put(channel, intensities.get(channel));
			}
		}
		return getAverage(toAverage, ignoreMissingValues);
	}

	public double getAverageFromNormalized(TurboIDExperimentType turboIDExperimentType, boolean ignoreMissingValues) {
		final TObjectDoubleHashMap<TurboID_Channel_Norm> intensities = getNormalizedIntensities(turboIDExperimentType);
		return getAverage(intensities, ignoreMissingValues);
	}

	public double getAverageFromOriginal(TurboIDExperimentType turboIDExperimentType, boolean ignoreMissingValues) {
		final TObjectDoubleHashMap<TurboID_Channel_Ori> intensities = getOriginalIntensities(turboIDExperimentType);
		return getAverage(intensities, ignoreMissingValues);
	}

	private double getAverage(TObjectDoubleHashMap<?> intensities, boolean ignoreMissingValues) {
		final TDoubleArrayList toAverage = new TDoubleArrayList();
		for (final double intensity : intensities.values()) {
			if (Double.isNaN(intensity)) {
				if (!ignoreMissingValues) {
					toAverage.add(0);
				}
			} else {
				toAverage.add(intensity);
			}
		}
		final double mean = Maths.mean(toAverage);

		return mean;
	}

	public String getGene() {
		return gene;
	}

	public boolean addOriginalIntensity(Double intensity, TurboID_Channel_Ori replicateID,
			TurboIDExperimentType turboIDExperiment) {
		if (Double.isNaN(intensity)) {
			return false;
		}
		if (originalIntensities.containsKey(turboIDExperiment)) {
			originalIntensities.get(turboIDExperiment).put(replicateID, intensity);
		} else {
			final TObjectDoubleHashMap<TurboID_Channel_Ori> intensitiesForThisProtein = new TObjectDoubleHashMap<TurboID_Channel_Ori>();
			intensitiesForThisProtein.put(replicateID, intensity);
			originalIntensities.put(turboIDExperiment, intensitiesForThisProtein);
		}
		return true;

	}

	public int getSpc(Replicate replicate) {
		return spcByReplicate.get(replicate);
	}

	public TObjectIntHashMap<Replicate> getSpcByReplicate() {
		return spcByReplicate;
	}

	public void setSpc(int spc, Replicate replicate) {
		spcByReplicate.put(replicate, spc);
	}

	public int getSumSPCAcrossReplicates() {
		if (sumSPCAcrossReplicates == -1) {
			sumSPCAcrossReplicates = 0;
			final Replicate[] values = Replicate.values();
			for (final Replicate replicate : values) {
				sumSPCAcrossReplicates += getSpc(replicate);
			}
		}
		return sumSPCAcrossReplicates;
	}

	public double getSumIntesitiesAcrossReplicates(TurboIDExperimentType expType, boolean useNormalizedIntensities) {
		if (useNormalizedIntensities) {
			final TObjectDoubleHashMap<TurboID_Channel_Norm> intensities = normalizedIntensities.get(expType);
			return new TDoubleArrayList(intensities.values()).sum();
		} else {
			final TObjectDoubleHashMap<TurboID_Channel_Ori> intensities = originalIntensities.get(expType);
			return new TDoubleArrayList(intensities.values()).sum();
		}
	}

	public double getStandardDeviationOfIntensities(TurboIDExperimentType expType, boolean useNormalizedIntensities) {
		if (useNormalizedIntensities) {
			final TObjectDoubleHashMap<TurboID_Channel_Norm> intensities = normalizedIntensities.get(expType);
			return Maths.stddev(intensities.values());
		} else {
			final TObjectDoubleHashMap<TurboID_Channel_Ori> intensities = originalIntensities.get(expType);
			return Maths.stddev(intensities.values());
		}
	}

	public double getAvgOfIntensities(TurboIDExperimentType expType, boolean useNormalizedIntensities) {
		if (useNormalizedIntensities) {
			final TObjectDoubleHashMap<TurboID_Channel_Norm> intensities = normalizedIntensities.get(expType);
			return Maths.mean(intensities.values());
		} else {
			final TObjectDoubleHashMap<TurboID_Channel_Ori> intensities = originalIntensities.get(expType);
			return Maths.mean(intensities.values());
		}
	}

	public double getStandardErrorOfMeasurementOfIntensities(TurboIDExperimentType expType,
			boolean useNormalizedIntensities) {
		if (useNormalizedIntensities) {
			final TObjectDoubleHashMap<TurboID_Channel_Norm> intensities = normalizedIntensities.get(expType);
			return Maths.sem(intensities.values());
		} else {
			final TObjectDoubleHashMap<TurboID_Channel_Ori> intensities = originalIntensities.get(expType);
			return Maths.sem(intensities.values());
		}
	}

	/**
	 * Gets the error measurement of the protein in one of the baits experiments
	 * as:<br>
	 * sqrt((sum(int)-N*avg(int))/(N-1))
	 * 
	 * @param expType
	 * @param useNormalizedIntensities
	 * @return
	 */
	public double getErrorMeasurement(TurboIDExperimentType expType, boolean useNormalizedIntensities) {
		double[] intensities;
		if (useNormalizedIntensities) {
			intensities = getNormalizedIntensities(expType).values();
		} else {
			intensities = getOriginalIntensities(expType).values();
		}
		if (getAcc().equals("Q64524")) {
			System.out.println(this);
		}
		final double sum = Maths.sum(intensities);
		final double avg = Maths.mean(intensities);
		final int n = intensities.length;
		double ret = Math.abs(sum - n * avg) / (1.0 * (n - 1));
		ret = Math.sqrt(ret);
		return ret;
	}

	public TObjectDoubleHashMap<TurboID_Channel_Ori> getPseudoSpecCountsWithOriginalIntensities(Replicate replicate) {
		if (!pseudoSpecCountsWithOriginalIntensitiesByReplicate.containsKey(replicate)) {
			final int spc = getSpc(replicate);
			final double sumIntesitiesAcrossChannels = getSumIntesitiesAcrossChannels(replicate, false);
			final TObjectDoubleHashMap<TurboID_Channel_Ori> ret = new TObjectDoubleHashMap<TurboID_Channel_Ori>();
			for (final TurboIDExperimentType expType : TurboIDExperimentType.values()) {

				final TObjectDoubleHashMap<TurboID_Channel_Ori> originalIntensities = getOriginalIntensities(expType);
				for (final TurboID_Channel_Ori channel : originalIntensities.keySet()) {
					if (channel.getReplicate() == replicate) {

						final double pseudoSPC = spc * originalIntensities.get(channel) / sumIntesitiesAcrossChannels;
						ret.put(channel, pseudoSPC);
					}
				}
			}
			pseudoSpecCountsWithOriginalIntensitiesByReplicate.put(replicate, ret);
			return ret;
		}
		return pseudoSpecCountsWithOriginalIntensitiesByReplicate.get(replicate);
	}

	public TObjectDoubleHashMap<TurboID_Channel_Norm> getPseudoSpecCountsWithNormalizedIntensities(
			Replicate replicate) {
		if (!pseudoSpecCountsWithNormalizedIntensitiesByReplicate.containsKey(replicate)) {

			final int spc = getSpc(replicate);
			final double sumIntesitiesAcrossChannels = getSumIntesitiesAcrossChannels(replicate, true);
			final TObjectDoubleHashMap<TurboID_Channel_Norm> ret = new TObjectDoubleHashMap<TurboID_Channel_Norm>();
			for (final TurboIDExperimentType expType : TurboIDExperimentType.values()) {
				// if (!expType.isBait()) {
				// continue;
				// }
				final TObjectDoubleHashMap<TurboID_Channel_Norm> normalizedIntensities = getNormalizedIntensities(
						expType);
				for (final TurboID_Channel_Norm channel : normalizedIntensities.keySet()) {
					if (channel.getReplicate() == replicate) {

						final double pseudoSPC = spc * normalizedIntensities.get(channel) / sumIntesitiesAcrossChannels;
						if (Double.compare(pseudoSPC, 0.0) == 0) {
							ret.put(channel, Double.NaN);
						} else {
							ret.put(channel, pseudoSPC);
						}
					}
				}
			}
			pseudoSpecCountsWithNormalizedIntensitiesByReplicate.put(replicate, ret);
			return ret;
		}
		return pseudoSpecCountsWithNormalizedIntensitiesByReplicate.get(replicate);
	}

	public TObjectDoubleHashMap<TurboID_Channel_Norm> getDistributedIntensitiesWithNormalizedIntensities() {
		if (distributedIntensitiesIntensities.isEmpty()) {

			for (final TurboIDExperimentType expType : TurboIDExperimentType.values()) {
				// if (!expType.isBait()) {
				// continue;
				// }

				final TObjectDoubleHashMap<TurboID_Channel_Norm> normalizedIntensities = getNormalizedIntensities(
						expType);
				for (final TurboID_Channel_Norm channel : normalizedIntensities.keySet()) {
					if (expType == TurboIDExperimentType.SUN1 && getGene().equals("Tbca")) {
						if (channel.name().equals("norm_Cy_A1")) {
							log.info("asdf");
						}
					}
					final double sumIntesitiesAcrossBaits = getSumIntesitiesAcrossChannels(channel.getReplicate(),
							true);
					final double d = normalizedIntensities.get(channel);
					// System.out.println(expType.name() + " " +
					// channel.name() + "\t" + d);
					if (Double.compare(d, 0.0) == 0) {
						distributedIntensitiesIntensities.put(channel, Double.NaN);
					} else {

						final double distributedIntensity = d / sumIntesitiesAcrossBaits;
						distributedIntensitiesIntensities.put(channel, distributedIntensity);
					}

				}
			}

		}
		return distributedIntensitiesIntensities;
	}

	private double getSumIntesitiesAcrossChannels(Replicate replicate, boolean useNormalizedIntensities) {
		if (useNormalizedIntensities) {
			double ret = 0.0;
			for (final TurboIDExperimentType expType : TurboIDExperimentType.values()) {

				final TObjectDoubleHashMap<TurboID_Channel_Norm> intensities = getNormalizedIntensities(expType);
				for (final TurboID_Channel_Norm channel : intensities.keySet()) {
					if (channel.getReplicate() == replicate) {
						final double intensity = intensities.get(channel);
						if (!Double.isNaN(intensity)) {
							ret += intensity;
						}
					}
				}
			}
			return ret;
		} else {
			double ret = 0.0;
			for (final TurboIDExperimentType expType : TurboIDExperimentType.values()) {
				final TObjectDoubleHashMap<TurboID_Channel_Ori> intensities = getOriginalIntensities(expType);
				for (final TurboID_Channel_Ori channel : intensities.keySet()) {
					if (channel.getReplicate() == replicate) {
						final double intensity = intensities.get(channel);
						if (!Double.isNaN(intensity)) {
							ret += intensity;
						}
					}
				}
			}
			return ret;
		}
	}

	private double getSumIntesitiesAcrossBaitChannels(boolean useNormalizedIntensities) {
		if (useNormalizedIntensities) {
			double ret = 0.0;
			for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaits()) {

				final TObjectDoubleHashMap<TurboID_Channel_Norm> intensities = getNormalizedIntensities(bait);
				for (final TurboID_Channel_Norm channel : intensities.keySet()) {

					final double intensity = intensities.get(channel);
					if (!Double.isNaN(intensity)) {
						ret += intensity;
					}

				}
			}
			return ret;
		} else {
			double ret = 0.0;
			for (final TurboIDExperimentType expType : TurboIDExperimentType.values()) {
				final TObjectDoubleHashMap<TurboID_Channel_Ori> intensities = getOriginalIntensities(expType);
				for (final TurboID_Channel_Ori channel : intensities.keySet()) {

					final double d = intensities.get(channel);
					if (!Double.isNaN(d)) {
						ret += d;
					}
				}
			}
			return ret;
		}
	}

	private double getMeanPseudoLog2SPC(TurboIDExperimentType bait, Replicate replicate) {
		final TObjectDoubleHashMap<TurboID_Channel_Norm> pseudoSpecCountsWithNormalizedIntensities = getPseudoSpecCountsWithNormalizedIntensities(
				replicate);
		final TDoubleArrayList toAverage = new TDoubleArrayList();
		for (final TurboID_Channel_Norm channel : pseudoSpecCountsWithNormalizedIntensities.keySet()) {
			if (channel.getExpType() == bait) {
				toAverage.add(pseudoSpecCountsWithNormalizedIntensities.get(channel));
			}
		}
		final double pseudoSPC = Maths.mean(toAverage);
		return Maths.log(pseudoSPC, 2);
	}

	private double getMeanPseudoLog2SPC(TurboIDExperimentType bait) {
		final TDoubleArrayList toAverage = new TDoubleArrayList();
		for (final Replicate replicate : Replicate.values()) {

			final TObjectDoubleHashMap<TurboID_Channel_Norm> pseudoSpecCountsWithNormalizedIntensities = getPseudoSpecCountsWithNormalizedIntensities(
					replicate);
			for (final TurboID_Channel_Norm channel : pseudoSpecCountsWithNormalizedIntensities.keySet()) {
				if (channel.getExpType() == bait) {
					toAverage.add(pseudoSpecCountsWithNormalizedIntensities.get(channel));
				}
			}
		}
		final double pseudoSPC = Maths.mean(toAverage);
		if (Double.compare(pseudoSPC, 0.0) == 0) {
			return Double.NaN;
		}
		return Maths.log(pseudoSPC, 2);
	}

	public boolean isSpecific(TurboIDExperimentType bait, Replicate replicate, double log2PseudoSPCThreshold) {

		final double pseudoLog2SPC = getMeanPseudoLog2SPC(bait, replicate);
		if (pseudoLog2SPC >= log2PseudoSPCThreshold) {
			return true;
		}
		return false;
	}

	public boolean isSpecificInBothExperiments(TurboIDExperimentType bait, double log2PseudoSPCThreshold_RepA,
			double log2PseudoSPCThreshold_RepB) {
		boolean ret = true;
		for (final Replicate replicate : Replicate.values()) {
			double threshold;
			if (replicate == Replicate.A) {
				threshold = log2PseudoSPCThreshold_RepA;
			} else {
				threshold = log2PseudoSPCThreshold_RepB;
			}
			ret = ret && isSpecific(bait, replicate, threshold);
		}
		return ret;
	}

	public double getLocalFDR(TurboIDExperimentType bait, Replicate replicate, Gaussian fit,
			List<WeightedObservedPoint> experimentalHistogram) {

		final double pseudoLog2SPC = getMeanPseudoLog2SPC(bait, replicate);
		final double fitY = fit.value(pseudoLog2SPC);
		double experimentalY = 0.0;
		for (int i = 0; i < experimentalHistogram.size(); i++) {
			final WeightedObservedPoint point = experimentalHistogram.get(i);
			if (point.getX() > pseudoLog2SPC) {
				if (i < experimentalHistogram.size() - 1) {
					experimentalY = Math
							.abs(experimentalHistogram.get(i + 1).getY() - experimentalHistogram.get(i).getY()) / 2;
				} else {
					experimentalY = experimentalHistogram.get(i).getY();
				}
				break;
			}
			// if we are in the last point and still pseudoLog2SPC is bigger,
			// then get the Y
			// value of the last point of the experimental histogram
			if (i == experimentalHistogram.size() - 1) {
				experimentalY = point.getY();
			}
		}
		final double fdr = fitY / experimentalY;
		return fdr;
	}

	public double getLocalFDR(TurboIDExperimentType bait, Gaussian fit,
			List<WeightedObservedPoint> experimentalHistogram) {

		final double pseudoLog2SPC = getMeanPseudoLog2SPC(bait);

		final double fitY = fit.value(pseudoLog2SPC);
		double experimentalY = 0.0;
		for (int i = 0; i < experimentalHistogram.size(); i++) {
			final WeightedObservedPoint point = experimentalHistogram.get(i);
			if (point.getX() > pseudoLog2SPC) {
				if (i < experimentalHistogram.size() - 1) {
					experimentalY = Math
							.abs(experimentalHistogram.get(i + 1).getY() - experimentalHistogram.get(i).getY()) / 2;
				} else {
					experimentalY = experimentalHistogram.get(i).getY();
				}
				break;
			}
			// if we are in the last point and still pseudoLog2SPC is bigger,
			// then get the Y
			// value of the last point of the experimental histogram
			if (i == experimentalHistogram.size() - 1) {
				experimentalY = point.getY();
			}
		}
		final double fdr = fitY / experimentalY;
		return fdr;
	}

	public double getAnovaPValueOverBaits(boolean applyLog, boolean distribSPC, boolean distribIntensity,
			TurboIDFraction fraction) {
		if (getGene().equalsIgnoreCase("cpe") && distribIntensity) {
			log.info(" ");
		}
		final List<double[]> inputAnova = new ArrayList<double[]>();
		inputAnova.clear();
		for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaits()) {
			final TDoubleArrayList intensities = new TDoubleArrayList();
			if (distribSPC) {
				for (final Replicate replicate : Replicate.values(fraction)) {

					final TObjectDoubleHashMap<TurboID_Channel_Norm> pseudoSpecCountsByChannel = getPseudoSpecCountsWithNormalizedIntensities(
							replicate);
					for (final TurboID_Channel_Norm channel : pseudoSpecCountsByChannel.keySet()) {
						if (channel.getExpType() == bait) {
							final double val = pseudoSpecCountsByChannel.get(channel);
							if (!Double.isNaN(val)) {
								intensities.add(val);
							}
						}
					}

				}

			} else if (distribIntensity) {

				final TObjectDoubleHashMap<TurboID_Channel_Norm> distributedIntensitiesByChannel = getDistributedIntensitiesWithNormalizedIntensities();
				for (final TurboID_Channel_Norm channel : distributedIntensitiesByChannel.keySet()) {
					if (channel.getExpType() == bait) {
						final double val = distributedIntensitiesByChannel.get(channel);
						if (!Double.isNaN(val)) {
							intensities.add(val);
						}
					}
				}

			} else {
				final TObjectDoubleHashMap<TurboID_Channel_Norm> normalizedIntensitiesByChannel = getNormalizedIntensities(
						bait);
				final double[] values = normalizedIntensitiesByChannel.values();
				for (final double d : values) {
					if (!Double.isNaN(d)) {
						intensities.add(d);
					}
				}

			}
			if (applyLog) {
				for (int i = 0; i < intensities.size(); i++) {
					intensities.set(i, Maths.log(intensities.get(i), 2));
				}
			}
			// 2 or more values are required per category
			if (intensities.size() > 1) {
				inputAnova.add(intensities.toArray());
			}
		}
		if (inputAnova.size() < 2) {
			return Double.NaN;
		}
		final OneWayAnova anova = new OneWayAnova();
		final double anovaPValue = anova.anovaPValue(inputAnova);

		return anovaPValue;
	}

	public boolean isAnovaValid(boolean applyLog, boolean distribSPC, boolean distribIntensity,
			TurboIDFraction fraction, int minPerGroup) {
		final List<double[]> inputAnova = new ArrayList<double[]>();
		inputAnova.clear();
		for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaits()) {
			final TDoubleArrayList intensities = new TDoubleArrayList();
			if (distribSPC) {
				for (final Replicate replicate : Replicate.values(fraction)) {

					final TObjectDoubleHashMap<TurboID_Channel_Norm> pseudoSpecCountsByChannel = getPseudoSpecCountsWithNormalizedIntensities(
							replicate);
					for (final TurboID_Channel_Norm channel : pseudoSpecCountsByChannel.keySet()) {
						if (channel.getExpType() == bait) {
							final double val = pseudoSpecCountsByChannel.get(channel);
							if (!Double.isNaN(val)) {
								intensities.add(val);
							}
						}
					}

				}

			} else if (distribIntensity) {

				final TObjectDoubleHashMap<TurboID_Channel_Norm> distributedIntensitiesByChannel = getDistributedIntensitiesWithNormalizedIntensities();
				for (final TurboID_Channel_Norm channel : distributedIntensitiesByChannel.keySet()) {
					if (channel.getExpType() == bait) {
						final double val = distributedIntensitiesByChannel.get(channel);
						if (!Double.isNaN(val)) {
							intensities.add(val);
						}
					}
				}

			} else {
				final TObjectDoubleHashMap<TurboID_Channel_Norm> normalizedIntensitiesByChannel = getNormalizedIntensities(
						bait);
				final double[] values = normalizedIntensitiesByChannel.values();
				for (final double d : values) {
					if (!Double.isNaN(d)) {
						intensities.add(d);
					}
				}

			}
			if (applyLog) {
				for (int i = 0; i < intensities.size(); i++) {
					intensities.set(i, Maths.log(intensities.get(i), 2));
				}
			}
			// 2 or more values are required per category
			if (intensities.size() > 1) {
				inputAnova.add(intensities.toArray());
			}
		}
		if (inputAnova.size() < 2) {
			return false;
		}

		for (final double[] ds : inputAnova) {
			if (ds.length < minPerGroup) {
				return false;
			}
		}

		return true;
	}

	public double getTTestPValueOverBaits(TurboIDExperimentType bait1, TurboIDExperimentType bait2, boolean applyLog,
			boolean distribSPC, boolean distribIntensity, TurboIDFraction fraction) {

		final List<double[]> inputTTest = new ArrayList<double[]>();
		inputTTest.clear();
		final List<TurboIDExperimentType> baits = new ArrayList<TurboIDExperimentType>();
		baits.add(bait1);
		baits.add(bait2);
		for (final TurboIDExperimentType bait : baits) {
			final TDoubleArrayList intensities = new TDoubleArrayList();
			if (distribSPC) {
				for (final Replicate replicate : Replicate.values(fraction)) {

					final TObjectDoubleHashMap<TurboID_Channel_Norm> pseudoSpecCountsByChannel = getPseudoSpecCountsWithNormalizedIntensities(
							replicate);
					for (final TurboID_Channel_Norm channel : pseudoSpecCountsByChannel.keySet()) {
						if (channel.getExpType() == bait) {
							final double val = pseudoSpecCountsByChannel.get(channel);
							if (!Double.isNaN(val)) {
								intensities.add(val);
							}
						}
					}

				}

			} else if (distribIntensity) {

				final TObjectDoubleHashMap<TurboID_Channel_Norm> distributedIntensitiesByChannel = getDistributedIntensitiesWithNormalizedIntensities();
				for (final TurboID_Channel_Norm channel : distributedIntensitiesByChannel.keySet()) {
					if (channel.getExpType() == bait) {
						final double val = distributedIntensitiesByChannel.get(channel);
						if (!Double.isNaN(val)) {
							intensities.add(val);
						}
					}
				}

			} else {
				final TObjectDoubleHashMap<TurboID_Channel_Norm> normalizedIntensitiesByChannel = getNormalizedIntensities(
						bait);
				final double[] values = normalizedIntensitiesByChannel.values();
				for (final double d : values) {
					if (!Double.isNaN(d)) {
						intensities.add(d);
					}
				}

			}
			if (applyLog) {
				for (int i = 0; i < intensities.size(); i++) {
					intensities.set(i, Maths.log(intensities.get(i), 2));
				}
			}
			// 2 or more values are required per category
			if (intensities.size() > 1) {
				inputTTest.add(intensities.toArray());
			}
		}
		if (inputTTest.size() < 2) {
			return Double.NaN;
		}
		final TTest ttest = new TTest();
		final double[] sample1 = inputTTest.get(0);
		final double[] sample2 = inputTTest.get(1);
		final double ttestPValue = ttest.tTest(sample1, sample2);
		final double ttestPValue2 = ttest.homoscedasticTTest(sample1, sample2);
//		log.info(ttestPValue + " - " + ttestPValue2);
		return ttestPValue2;
	}

	/**
	 * Same as getAnovaPValueOverBaits, but including in the test the channel
	 * TurboID
	 * 
	 * @param applyLog
	 * @param distribSPC
	 * @param distribIntensity
	 * @param fraction
	 * @return
	 */
	public double getAnovaPValueOverBaitsAndTurboIDOnly(boolean applyLog, boolean distribSPC, boolean distribIntensity,
			TurboIDFraction fraction) {

		final List<double[]> inputAnova = new ArrayList<double[]>();
		inputAnova.clear();
		for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaitsAndTBID()) {
			final TDoubleArrayList intensities = new TDoubleArrayList();
			if (distribSPC) {
				for (final Replicate replicate : Replicate.values(fraction)) {

					final TObjectDoubleHashMap<TurboID_Channel_Norm> pseudoSpecCountsByChannel = getPseudoSpecCountsWithNormalizedIntensities(
							replicate);
					for (final TurboID_Channel_Norm channel : pseudoSpecCountsByChannel.keySet()) {
						if (channel.getExpType() == bait) {
							final double val = pseudoSpecCountsByChannel.get(channel);
							if (!Double.isNaN(val)) {
								intensities.add(val);
							}
						}
					}

				}

			} else if (distribIntensity) {

				final TObjectDoubleHashMap<TurboID_Channel_Norm> distributedIntensitiesByChannel = getDistributedIntensitiesWithNormalizedIntensities();
				for (final TurboID_Channel_Norm channel : distributedIntensitiesByChannel.keySet()) {
					if (channel.getExpType() == bait) {
						final double val = distributedIntensitiesByChannel.get(channel);
						if (!Double.isNaN(val)) {
							intensities.add(val);
						}
					}
				}

			} else {
				final TObjectDoubleHashMap<TurboID_Channel_Norm> normalizedIntensitiesByChannel = getNormalizedIntensities(
						bait);
				final double[] values = normalizedIntensitiesByChannel.values();
				for (final double d : values) {
					if (!Double.isNaN(d)) {
						intensities.add(d);
					}
				}

			}
			if (applyLog) {
				for (int i = 0; i < intensities.size(); i++) {
					intensities.set(i, Maths.log(intensities.get(i), 2));
				}
			}
			// 2 or more values are required per category
			if (intensities.size() > 1) {
				inputAnova.add(intensities.toArray());
			}
		}
		if (inputAnova.size() < 2) {
			return Double.NaN;
		}
		final OneWayAnova anova = new OneWayAnova();
		final double anovaPValue = anova.anovaPValue(inputAnova);

		return anovaPValue;
	}

}
