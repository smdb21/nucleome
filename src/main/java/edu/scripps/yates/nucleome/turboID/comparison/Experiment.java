package edu.scripps.yates.nucleome.turboID.comparison;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.ml.distance.EuclideanDistance;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

import edu.scripps.yates.nucleome.turboID.TurboIDExperimentType;
import edu.scripps.yates.utilities.maths.Maths;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;

public class Experiment {
	private final TIntObjectHashMap<List<ProteinToCompare>> proteinsByCluster = new TIntObjectHashMap<List<ProteinToCompare>>();
	private final TIntObjectHashMap<TObjectDoubleMap<TurboIDExperimentType>> avgQuantsByCluster = new TIntObjectHashMap<TObjectDoubleMap<TurboIDExperimentType>>();
	private final Map<String, ProteinToCompare> proteinsByAcc = new THashMap<String, ProteinToCompare>();
	private final List<ProteinToCompare> proteins = new ArrayList<ProteinToCompare>();
	private int totalSPC;
	private Set<String> totalProteins;
	private int numCompleteProteins;
	private Set<String> clusteredProteins;

	public Experiment() {

	}

	public void addProtein(ProteinToCompare protein) {
		if (protein.getClusterNum() != null) {
			if (!proteinsByCluster.containsKey(protein.getClusterNum())) {
				proteinsByCluster.put(protein.getClusterNum(), new ArrayList<ProteinToCompare>());
			}
			proteinsByCluster.get(protein.getClusterNum()).add(protein);
		}
		proteinsByAcc.put(protein.getAcc(), protein);
		proteins.add(protein);
	}

	public int getNumSize() {
		return proteins.size();
	}

	public List<ProteinToCompare> getProteinsByClusterNumber(int cluster) {
		return proteinsByCluster.get(cluster);
	}

	public TObjectDoubleMap<TurboIDExperimentType> getAverageQuantsByCluster(int cluster) {
		if (!avgQuantsByCluster.containsKey(cluster)) {
			final TObjectDoubleMap<TurboIDExperimentType> ret = new TObjectDoubleHashMap<TurboIDExperimentType>();
			final List<ProteinToCompare> proteinList = proteinsByCluster.get(cluster);
			for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaits()) {
				final TDoubleList toAverage = new TDoubleArrayList();
				for (final ProteinToCompare protein : proteinList) {
					final double quant = protein.getQuant(bait);
					if (!Double.isNaN(quant)) {
						toAverage.add(quant);
					}
				}
				final double avg = Maths.mean(toAverage);
				ret.put(bait, avg);
			}
			avgQuantsByCluster.put(cluster, ret);
		}
		return avgQuantsByCluster.get(cluster);
	}

	public int getNumClusters() {
		return proteinsByCluster.size();
	}

	private double[] getArray(TObjectDoubleMap<TurboIDExperimentType> averageQuantsByCluster) {
		final double[] ret = new double[averageQuantsByCluster.size()];
		int index = 0;
		for (final TurboIDExperimentType bait : TurboIDExperimentType.getBaits()) {
			if (averageQuantsByCluster.containsKey(bait)) {
				ret[index++] = averageQuantsByCluster.get(bait);
			}
		}
		return ret;
	}

	public int getCloserClusterByEuclideanDistance(TObjectDoubleMap<TurboIDExperimentType> averageQuantsByCluster) {

		final double[] a = getArray(averageQuantsByCluster);

		Double min = Double.MAX_VALUE;
		int clusterNum = -1;
		for (final int cluster2 : proteinsByCluster.keys()) {
			final double[] b = getArray(getAverageQuantsByCluster(cluster2));
			final EuclideanDistance distance = new EuclideanDistance();
			final double dist = distance.compute(a, b);
			if (dist < min) {
				min = dist;
				clusterNum = cluster2;
			}
		}
		return clusterNum;
	}

	public int getCloserClusterByPearsonCorrelation(TObjectDoubleMap<TurboIDExperimentType> averageQuantsByCluster) {

		final double[] a = getArray(averageQuantsByCluster);

		Double max = -Double.MAX_VALUE;
		int clusterNum = -1;
		for (final int cluster2 : proteinsByCluster.keys()) {
			final double[] b = getArray(getAverageQuantsByCluster(cluster2));
			final PearsonsCorrelation correlation = new PearsonsCorrelation();
			final double dist = correlation.correlation(a, b);
			if (dist > max) {
				max = dist;
				clusterNum = cluster2;
			}
		}
		return clusterNum;
	}

	public List<ProteinToCompare> getProteins() {
		return proteins;
	}

	public Map<String, ProteinToCompare> getProteinsByAcc() {
		return proteinsByAcc;
	}

	public void setTotalSPC(int totalSPCInExp) {
		totalSPC = totalSPCInExp;

	}

	public int getTotalSPC() {
		return totalSPC;
	}

	public Set<String> getTotalProteins() {
		return totalProteins;
	}

	public void setTotalProteins(Set<String> proteins2) {
		totalProteins = proteins2;
	}

	public void setNumCompleteProteins(int numValid) {
		numCompleteProteins = numValid;

	}

	public int getNumCompleteProteins() {
		return numCompleteProteins;
	}

	public void setClusteredProteins(Set<String> clusteredProteins) {
		this.clusteredProteins = clusteredProteins;
	}

	public Set<String> getClusteredProteins() {
		return clusteredProteins;
	}

}