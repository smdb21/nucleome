package edu.scripps.yates.nucleome.model;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.nucleome.Constants;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.proteomicsmodel.Protein;
import edu.scripps.yates.utilities.remote.RemoteSSHFileReference;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;

public class Experiment {
	private final String experimentName;
	private final CellType cellType;
	private final Map<Integer, Replicate> replicates = new HashMap<Integer, Replicate>();
	private Set<String> proteinAccs;
	private final Wash wash;

	public Experiment(String name, Wash wash, CellType cellType) {
		experimentName = name;
		this.cellType = cellType;
		this.wash = wash;
	}

	public Experiment(String name, CellType cellType) {
		this(name, null, cellType);
	}

	public void addReplicate(int num, CellType cellType, CellCompartment cellCompartment, File file)
			throws IOException {
		addReplicate(num, null, cellType, cellCompartment, file);
	}

	public void addReplicate(int num, Wash wash, CellType cellType, CellCompartment cellCompartment, File file)
			throws IOException {
		addReplicate(num, wash, cellType, cellCompartment, new RemoteSSHFileReference(file));
	}

	public void addReplicate(int replicateNum, Wash wash, CellType cellType, CellCompartment cellCompartment,
			RemoteSSHFileReference remote) throws IOException {
		if (cellType != this.cellType) {
			throw new IllegalArgumentException("Different cell types");
		}
		Replicate replicate = null;
		if (replicates.containsKey(replicateNum)) {
			replicate = replicates.get(replicateNum);
		} else {
			replicate = new Replicate(experimentName, replicateNum, wash, cellType);
			replicates.put(replicateNum, replicate);
		}
		replicate.setFraction(cellCompartment, wash, remote);
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

		final TIntArrayList values = new TIntArrayList();
		for (final Replicate replicate : replicates.values()) {
			final Fractionation fractionation = replicate.getFractionation(cellCompartment);
			if (fractionation != null) {
				values.add(fractionation.getSpectralCount(proteinAcc, skipFilters));
			}
		}
		if (!values.isEmpty()) {
			return Maths.mean(values);
		}
		return 0;
	}

	public double getAvgPeptideCount(String proteinAcc, CellCompartment cellCompartment, boolean skipFilters)
			throws IOException {

		final TIntArrayList values = new TIntArrayList();
		for (final Replicate replicate : replicates.values()) {
			final Fractionation fractionation = replicate.getFractionation(cellCompartment);
			if (fractionation != null) {
				values.add(fractionation.getPeptideCount(proteinAcc, skipFilters));
			}
		}
		if (!values.isEmpty()) {
			return Maths.mean(values);
		}
		return 0.0;
	}

	public double getAvgSpectralCount(Collection<String> proteinAccessions, CellCompartment cellCompartment,
			boolean skipFilters) throws IOException {

		final TIntArrayList values = new TIntArrayList();
		for (final Replicate replicate : replicates.values()) {
			final Fractionation fractionation = replicate.getFractionation(cellCompartment);
			if (fractionation != null) {
				values.add(fractionation.getSpectralCount(proteinAccessions, skipFilters));
			}
		}
		if (!values.isEmpty()) {
			return Maths.mean(values);
		}
		return 0.0;
	}

	public double getAvgNSAF(Collection<String> proteinAccessions, CellCompartment cellCompartment, boolean skipFilters)
			throws IOException {
		final TDoubleArrayList values = new TDoubleArrayList();
		for (final Replicate replicate : replicates.values()) {
			final Fractionation fractionation = replicate.getFractionation(cellCompartment);
			if (fractionation != null) {
				values.add(fractionation.getAverageNSAF(proteinAccessions, skipFilters));
			}
		}
		if (!values.isEmpty()) {
			return Maths.mean(values);
		}
		return 0.0;
	}

	public double getAvgPeptideCount(Collection<String> proteinAccessions, CellCompartment cellCompartment,
			boolean skipFilters) throws IOException {

		final TIntArrayList values = new TIntArrayList();
		for (final Replicate replicate : replicates.values()) {
			final Fractionation fractionation = replicate.getFractionation(cellCompartment);
			if (fractionation != null) {
				values.add(fractionation.getPeptideCount(proteinAccessions, skipFilters));
			}
		}
		if (!values.isEmpty()) {
			return Maths.mean(values);
		}
		return 0;
	}

	public double getSumSPC(String proteinAcc, CellCompartment cellCompartment, boolean skipFilters)
			throws IOException {
		double ret = 0.0;
		for (final Replicate replicate : replicates.values()) {

			final Fractionation fractionation = replicate.getFractionation(cellCompartment);
			if (fractionation != null) {
				ret += fractionation.getSpectralCount(proteinAcc, skipFilters);
			}
		}
		return ret;
	}

	public int getSumSPC(Collection<String> proteinAccessions, CellCompartment cellCompartment, boolean skipFilters)
			throws IOException {
		int ret = 0;
		for (final Replicate replicate : replicates.values()) {

			final Fractionation fractionation = replicate.getFractionation(cellCompartment);
			if (fractionation != null) {
				ret += fractionation.getSpectralCount(proteinAccessions, skipFilters);
			}
		}
		return ret;
	}

	public double getSumNSAF(Collection<String> proteinAccessions, CellCompartment cellCompartment, boolean skipFilters)
			throws IOException {
		double ret = 0;
		for (final Replicate replicate : replicates.values()) {

			final Fractionation fractionation = replicate.getFractionation(cellCompartment);
			if (fractionation != null) {
				ret += fractionation.getSumNSAF(proteinAccessions, skipFilters);
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

	public double getSPCRatio(Collection<String> proteinAccessions, CellCompartment cellCompartment1,
			CellCompartment cellCompartment2, boolean skipFilters) throws IOException {

		double spc1 = getSumSPC(proteinAccessions, cellCompartment1, skipFilters);
		if (!skipFilters) {
			if (getAvgSpectralCount(proteinAccessions, cellCompartment1, true) < Constants.MIN_AVG_SPC) {
				spc1 = 0.0;
			}
		}
		double spc2 = getSumSPC(proteinAccessions, cellCompartment2, skipFilters);
		if (!skipFilters) {
			if (getAvgSpectralCount(proteinAccessions, cellCompartment2, true) < Constants.MIN_AVG_SPC) {
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
		final double nonLogValue = spc1 / spc2;
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

		final Set<String> ret = new HashSet<String>();
		for (final Replicate replicate : replicates.values()) {
			final Fractionation fractionation = replicate.getFractionation(cellCompartment);
			if (fractionation != null) {
				ret.addAll(fractionation.getProteinAccs());
			}
		}

		return ret;
	}

	public Set<Protein> getProteins(CellCompartment cellCompartment) throws IOException {

		final Set<Protein> ret = new HashSet<Protein>();
		for (final Replicate replicate : replicates.values()) {
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
			for (final CellCompartment cellCompartment : CellCompartment.values()) {
				proteinAccs.addAll(getProteinAccs(cellCompartment));
			}
		}
		return proteinAccs;
	}

	public List<Replicate> getReplicates() {
		final List<Replicate> ret = new ArrayList<Replicate>();
		for (int i = 1; i <= replicates.size(); i++) {
			ret.add(replicates.get(i));
		}
		return ret;
	}

	public Map<Integer, Replicate> getReplicateMap() {
		return replicates;
	}

	public String printHeader(boolean writeMSRunPresence) {
		final StringBuilder sb = new StringBuilder();
		for (final Replicate replicate : getSortedReplicates()) {
			sb.append(replicate.printHeader());
		}
		if (writeMSRunPresence) {
			sb.append(getName() + "_MSRunsPresence\t");
		}
		return sb.toString();
	}

	public String printColumnsForProtein(String rawAcc, boolean peptideCount) throws IOException {
		final StringBuilder sb = new StringBuilder();
		for (final Replicate replicate : getSortedReplicates()) {
			sb.append(replicate.printColumnsForProtein(rawAcc, peptideCount));
		}
		return sb.toString();
	}

	public String printColumnsForProteinGroup(Collection<String> proteinAccessions, DataType dataType)
			throws IOException {
		final StringBuilder sb = new StringBuilder();
		for (final Replicate replicate : getSortedReplicates()) {
			sb.append(replicate.printColumnsForProteinGroup(proteinAccessions, dataType));
		}

		return sb.toString();
	}

	public List<Replicate> getSortedReplicates() {
		final List<Replicate> list = new ArrayList<Replicate>();
		for (final Replicate replicate : replicates.values()) {
			list.add(replicate);
		}
		Collections.sort(list, new Comparator<Replicate>() {

			@Override
			public int compare(Replicate o1, Replicate o2) {
				return Integer.compare(o1.getReplicateNum(), o2.getReplicateNum());
			}

		});
		return list;
	}

	public Set<Protein> getProteins() throws IOException {
		final Set<Protein> ret = new HashSet<Protein>();

		for (final CellCompartment cellCompartment : CellCompartment.values()) {
			ret.addAll(getProteins(cellCompartment));
		}

		return ret;
	}

	public int getNumReplicatesWithSPCGreaterThan(Collection<String> proteinAccessions, CellCompartment cellCompartment,
			int spcThreshold) throws IOException {
		int ret = 0;
		final List<Replicate> sortedReplicates = getSortedReplicates();
		for (final Replicate replicate : sortedReplicates) {
			final Fractionation fractionation = replicate.getFractionation(cellCompartment);
			if (fractionation != null) {
				final int spc = fractionation.getSpectralCount(proteinAccessions, true);
				if (spc > spcThreshold) {
					ret++;
				}
			}
		}
		return ret;
	}

	public Wash getWash() {
		return wash;
	}

	@Override
	public String toString() {
		return getName();
	}
}
