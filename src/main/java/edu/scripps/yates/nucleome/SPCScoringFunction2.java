package edu.scripps.yates.nucleome;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;

import edu.scripps.yates.nucleome.model.CellCompartment;
import edu.scripps.yates.nucleome.model.CellType;
import edu.scripps.yates.nucleome.model.Experiment;
import edu.scripps.yates.utilities.grouping.ProteinGroup;

public class SPCScoringFunction2 extends ScoringFunction {
	private final static Logger log = Logger.getLogger(SPCScoringFunction2.class);
	private final _4DNucleomeAnalyzer analyzer;

	public SPCScoringFunction2(_4DNucleomeAnalyzer analyzer) {
		this.analyzer = analyzer;
	}

	@Override
	public double getScore(ProteinGroup proteinGroup, CellType celltype) throws IOException {
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
				double spcToStudy = experiment.getAvgSpectralCount(proteinGroup, Constants.cellCompartmentToStudy,
						true);
				double spcDenominator = experiment.getAvgSpectralCount(proteinGroup, denominatorCellCompartment, true);
				if (spcToStudy > 0 && spcDenominator > 0) {
					ratioInThisExperiment = true;
				}
				if (spcToStudy >= spcDenominator) {
					score += 1;
					if (spcToStudy >= 3 * spcDenominator && spcToStudy > 3) {
						score += 1;
					}
				} else {
					if (Constants.includeNegativeScoring) {
						score -= 1;
						if (3 * spcToStudy <= spcDenominator && spcDenominator > 3) {
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
