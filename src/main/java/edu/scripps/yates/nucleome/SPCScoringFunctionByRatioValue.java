package edu.scripps.yates.nucleome;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.apache.log4j.Logger;

import edu.scripps.yates.nucleome.model.CellCompartment;
import edu.scripps.yates.nucleome.model.CellType;
import edu.scripps.yates.nucleome.model.Experiment;

public class SPCScoringFunctionByRatioValue extends ScoringFunction {
	private final static Logger log = Logger.getLogger(SPCScoringFunctionByRatioValue.class);
	private final _4DNucleomeAnalyzer analyzer;

	public SPCScoringFunctionByRatioValue(_4DNucleomeAnalyzer analyzer) {
		this.analyzer = analyzer;
	}

	@Override
	public double getScore(Collection<String> proteinAccessions, CellType celltype) throws IOException {
		int numExperiments = 0;
		double score = 0.0;
		List<Experiment> experimentList = new ArrayList<Experiment>();
		if (celltype == null || celltype == CellType.A) {
			experimentList.addAll(analyzer.getExperimentsA());
		}
		if (celltype == null || celltype == CellType.M) {
			experimentList.addAll(analyzer.getExperimentsM());
		}
		if (celltype == null || celltype == CellType.U) {
			experimentList.addAll(analyzer.getExperimentsU());
		}

		// log.info(experimentList.size() + " experiments");
		for (Experiment experiment : experimentList) {
			boolean ratioInThisExperiment = false;

			for (CellCompartment denominatorCellCompartment : CellCompartment.values()) {
				if (denominatorCellCompartment == Constants.cellCompartmentToStudy) {
					continue;
				}

				// NE/N
				double spcRatio = experiment.getSPCRatio(proteinAccessions, Constants.cellCompartmentToStudy,
						denominatorCellCompartment, false);
				if (!Double.isNaN(spcRatio)) {
					ratioInThisExperiment = true;
					double log2Ratio = getLog2Ratio(spcRatio);
					if (log2Ratio > 10) {
						score += 1;
					}
					if (log2Ratio > 3) {
						score += 1;
					}
					if (log2Ratio > 0) {
						score += 1;
					}
					if (Constants.includeNegativeScoring) {
						if (log2Ratio < 0) {
							score -= 1;
						}
						if (log2Ratio < -3) {
							score -= 1;
						}
						if (log2Ratio < -10) {
							score -= 1;
						}
					}
				}

			}
			if (ratioInThisExperiment) {
				numExperiments++;
			}

		}
		// if we are looking in a particular cellType, divide
		score = score / numExperiments;
		return score;
	}
}
