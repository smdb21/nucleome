package edu.scripps.yates.nucleome.model;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.nucleome.Constants;
import edu.scripps.yates.utilities.grouping.ProteinGroup;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.proteomicsmodel.Protein;
import edu.scripps.yates.utilities.remote.RemoteSSHFileReference;

public class Experiment {
	private final String experimentName;
	private final CellType cellType;
	private final Map<Integer, Replicate> replicates = new HashMap<Integer, Replicate>();
	private Set<String> proteinAccs;

	public Experiment(String name, CellType cellType) {
		experimentName = name;
		this.cellType = cellType;
	}

	public void addReplicate(int num, CellType cellType, CellCompartment cellCompartment, File file)
			throws IOException {
		addReplicate(num, cellType, cellCompartment, new RemoteSSHFileReference(file));
	}

	public void addReplicate(int replicateNum, CellType cellType, CellCompartment cellCompartment,
			RemoteSSHFileReference remote) throws IOException {
		if (cellType != this.cellType) {
			throw new IllegalArgumentException("Different cell types");
		}
		Replicate replicate = null;
		if (replicates.containsKey(replicateNum)) {
			replicate = replicates.get(replicateNum);
		} else {
			replicate = new Replicate(experimentName, replicateNum, cellType);
			replicates.put(replicateNum, replicate);
		}
		replicate.setFraction(cellCompartment, remote);
	}

	/**
	 * @return the name
	 */
	public String getName() {
		return experimentName;
	}

	/**
	 * @return the cellType
	 */
	public CellType getCellType() {
		return cellType;
	}

	public double getAvgSpectralCount(String proteinAcc, CellCompartment cellCompartment, boolean skipFilters)
			throws IOException {
		List<Integer> values = new ArrayList<Integer>();
		for (Replicate replicate : replicates.values()) {
			final Fractionation fractionation = replicate.getFractionation(cellCompartment);
			if (fractionation != null) {
				values.add(fractionation.getSpectralCount(proteinAcc, skipFilters));
			}
		}
		if (!values.isEmpty()) {
			return Maths.mean(values.toArray(new Integer[0]));
		}
		return 0;
	}

	public double getAvgPeptideCount(String proteinAcc, CellCompartment cellCompartment, boolean skipFilters)
			throws IOException {
		List<Integer> values = new ArrayList<Integer>();
		for (Replicate replicate : replicates.values()) {
			final Fractionation fractionation = replicate.getFractionation(cellCompartment);
			if (fractionation != null) {
				values.add(fractionation.getPeptideCount(proteinAcc, skipFilters));
			}
		}
		if (!values.isEmpty()) {
			return Maths.mean(values.toArray(new Integer[0]));
		}
		return 0.0;
	}

	public double getAvgSpectralCount(ProteinGroup proteinGroup, CellCompartment cellCompartment, boolean skipFilters)
			throws IOException {
		List<Integer> values = new ArrayList<Integer>();
		for (Replicate replicate : replicates.values()) {
			final Fractionation fractionation = replicate.getFractionation(cellCompartment);
			if (fractionation != null) {
				values.add(fractionation.getSpectralCount(proteinGroup, skipFilters));
			}
		}
		if (!values.isEmpty()) {
			return Maths.mean(values.toArray(new Integer[0]));
		}
		return 0.0;
	}

	public double getAvgPeptideCount(ProteinGroup proteinGroup, CellCompartment cellCompartment, boolean skipFilters)
			throws IOException {
		List<Integer> values = new ArrayList<Integer>();
		for (Replicate replicate : replicates.values()) {
			final Fractionation fractionation = replicate.getFractionation(cellCompartment);
			if (fractionation != null) {
				values.add(fractionation.getPeptideCount(proteinGroup, skipFilters));
			}
		}
		if (!values.isEmpty()) {
			return Maths.mean(values.toArray(new Integer[0]));
		}
		return 0;
	}

	public double getSumSPC(String proteinAcc, CellCompartment cellCompartment, boolean skipFilters)
			throws IOException {
		double ret = 0.0;
		for (Replicate replicate : replicates.values()) {

			final Fractionation fractionation = replicate.getFractionation(cellCompartment);
			if (fractionation != null) {
				ret += fractionation.getSpectralCount(proteinAcc, skipFilters);
			}
		}
		return ret;
	}

	public int getSumSPC(ProteinGroup proteinGroup, CellCompartment cellCompartment, boolean skipFilters)
			throws IOException {
		int ret = 0;
		for (Replicate replicate : replicates.values()) {

			final Fractionation fractionation = replicate.getFractionation(cellCompartment);
			if (fractionation != null) {
				ret += fractionation.getSpectralCount(proteinGroup, skipFilters);
			}
		}
		return ret;
	}

	public double getSPCRatio(String proteinAcc, CellCompartment cellCompartment1, CellCompartment cellCompartment2,
			boolean skipFilters) throws IOException {
		double spc1 = getSumSPC(proteinAcc, cellCompartment1, skipFilters);
		if (!skipFilters) {
			if (getAvgSpectralCount(proteinAcc, cellCompartment1, true) < Constants.MIN_AVG_SPC) {
				spc1 = 0.0;
			}
		}
		double spc2 = getSumSPC(proteinAcc, cellCompartment2, skipFilters);
		if (!skipFilters) {
			if (getAvgSpectralCount(proteinAcc, cellCompartment2, true) < Constants.MIN_AVG_SPC) {
				spc2 = 0.0;
			}
		}
		return getNonLog2Ratio(spc1, spc2);
	}

	public double getSPCRatio(ProteinGroup proteinGroup, CellCompartment cellCompartment1,
			CellCompartment cellCompartment2, boolean skipFilters) throws IOException {

		double spc1 = getSumSPC(proteinGroup, cellCompartment1, skipFilters);
		if (!skipFilters) {
			if (getAvgSpectralCount(proteinGroup, cellCompartment1, true) < Constants.MIN_AVG_SPC) {
				spc1 = 0.0;
			}
		}
		double spc2 = getSumSPC(proteinGroup, cellCompartment2, skipFilters);
		if (!skipFilters) {
			if (getAvgSpectralCount(proteinGroup, cellCompartment2, true) < Constants.MIN_AVG_SPC) {
				spc2 = 0.0;
			}
		}
		return getNonLog2Ratio(spc1, spc2);
	}

	private double getLog2Ratio(double spc1, double spc2) {
		if (spc1 == 0.0 && spc2 == 0.0) {
			return Double.NaN;
		}
		if (spc1 == 0.0) {
			return Double.POSITIVE_INFINITY;
		}
		if (spc2 == 0.0) {
			return Double.NEGATIVE_INFINITY;
		}
		double nonLogValue = spc1 / spc2;
		return Math.log(nonLogValue) / Math.log(2);
	}

	private double getNonLog2Ratio(double spc1, double spc2) {
		if (spc1 == 0.0 && spc2 == 0.0) {
			return Double.NaN;
		}

		if (spc2 == 0.0) {
			return Double.POSITIVE_INFINITY;
		}
		return spc1 / spc2;
	}

	public Set<String> getProteinAccs(CellCompartment cellCompartment) throws IOException {

		Set<String> ret = new HashSet<String>();
		for (Replicate replicate : replicates.values()) {
			final Fractionation fractionation = replicate.getFractionation(cellCompartment);
			if (fractionation != null) {
				ret.addAll(fractionation.getProteinAccs());
			}
		}

		return ret;
	}

	public Set<Protein> getProteins(CellCompartment cellCompartment) throws IOException {

		Set<Protein> ret = new HashSet<Protein>();
		for (Replicate replicate : replicates.values()) {
			final Fractionation fractionation = replicate.getFractionation(cellCompartment);
			if (fractionation != null) {
				ret.addAll(fractionation.getProteins());
			}
		}

		return ret;
	}

	public Set<String> getProteinAccs() throws IOException {
		if (proteinAccs == null) {
			proteinAccs = new HashSet<String>();
			for (CellCompartment cellCompartment : CellCompartment.values()) {
				proteinAccs.addAll(getProteinAccs(cellCompartment));
			}
		}
		return proteinAccs;
	}

	public List<Replicate> getReplicates() {
		List<Replicate> ret = new ArrayList<Replicate>();
		for (int i = 1; i < replicates.size(); i++) {
			ret.add(replicates.get(i));
		}
		return ret;
	}

	public Map<Integer, Replicate> getReplicateMap() {
		return replicates;
	}

	public String printHeader() {
		StringBuilder sb = new StringBuilder();
		for (Replicate replicate : getSortedReplicates()) {
			sb.append(replicate.printHeader());
		}
		return sb.toString();
	}

	public String printColumnsForProtein(String rawAcc, boolean peptideCount) throws IOException {
		StringBuilder sb = new StringBuilder();
		for (Replicate replicate : getSortedReplicates()) {
			sb.append(replicate.printColumnsForProtein(rawAcc, peptideCount));
		}
		return sb.toString();
	}

	public String printColumnsForProteinGroup(ProteinGroup proteinGroup, boolean peptideCount) throws IOException {
		StringBuilder sb = new StringBuilder();
		for (Replicate replicate : getSortedReplicates()) {
			sb.append(replicate.printColumnsForProteinGroup(proteinGroup, peptideCount));
		}
		return sb.toString();
	}

	private List<Replicate> getSortedReplicates() {
		List<Replicate> list = new ArrayList<Replicate>();
		for (Replicate replicate : replicates.values()) {
			list.add(replicate);
		}
		Collections.sort(list, new Comparator<Replicate>() {

			public int compare(Replicate o1, Replicate o2) {
				return Integer.compare(o1.getReplicateNum(), o2.getReplicateNum());
			}

		});
		return list;
	}

	public Set<Protein> getProteins() throws IOException {
		Set<Protein> ret = new HashSet<Protein>();

		for (CellCompartment cellCompartment : CellCompartment.values()) {
			ret.addAll(getProteins(cellCompartment));
		}

		return ret;
	}
}
