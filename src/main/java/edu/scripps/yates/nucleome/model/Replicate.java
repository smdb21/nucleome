package edu.scripps.yates.nucleome.model;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import edu.scripps.yates.utilities.grouping.ProteinGroup;
import edu.scripps.yates.utilities.remote.RemoteSSHFileReference;

/**
 * This would be like U1 (rep1), having N, Ne and C
 *
 * @author Salva
 *
 */
public class Replicate {
	private final int replicateNum;
	private final CellType cellType;
	private final String experimentName;
	private final Map<CellCompartment, Fractionation> fractions = new HashMap<CellCompartment, Fractionation>();

	public Replicate(String experimentName, int replicateNum, CellType cellType) {
		this.replicateNum = replicateNum;
		this.cellType = cellType;
		this.experimentName = experimentName;
	}

	public void setFraction(CellCompartment cellCompartment, RemoteSSHFileReference remote) throws IOException {
		Fractionation fraction = new Fractionation(experimentName, replicateNum, cellCompartment, remote, cellType);
		fractions.put(cellCompartment, fraction);
	}

	/**
	 * @return the num
	 */
	public int getReplicateNum() {
		return replicateNum;
	}

	public Fractionation getFractionation(CellCompartment cellCompartment) {
		return fractions.get(cellCompartment);
	}

	public String printColumnsForProtein(String rawAcc, boolean peptideCount) throws IOException {
		StringBuilder sb = new StringBuilder();
		for (CellCompartment cellCompartment : CellCompartment.values()) {
			final Fractionation fractionation = getFractionation(cellCompartment);
			if (fractionation != null) {
				if (peptideCount) {
					sb.append(fractionation.getPeptideCount(rawAcc, true));
				} else {
					sb.append(fractionation.getSpectralCount(rawAcc, true));
				}
				sb.append("\t");
			}

		}
		return sb.toString();
	}

	public String printColumnsForProteinGroup(ProteinGroup proteinGroup, boolean peptideCount) throws IOException {
		StringBuilder sb = new StringBuilder();
		for (CellCompartment cellCompartment : CellCompartment.values()) {
			final Fractionation fractionation = getFractionation(cellCompartment);
			if (fractionation != null) {
				if (peptideCount) {
					sb.append(fractionation.getPeptideCount(proteinGroup, true));
				} else {
					sb.append(fractionation.getSpectralCount(proteinGroup, true));
				}
				sb.append("\t");
			}

		}
		return sb.toString();
	}

	public String printHeader() {
		StringBuilder sb = new StringBuilder();
		for (CellCompartment cellCompartment : CellCompartment.values()) {
			final Fractionation fractionation = getFractionation(cellCompartment);
			if (fractionation != null) {
				sb.append(fractionation.getName());
				sb.append("\t");
			}

		}
		return sb.toString();
	}

}
