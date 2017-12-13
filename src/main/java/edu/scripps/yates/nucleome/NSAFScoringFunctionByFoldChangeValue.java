package edu.scripps.yates.nucleome;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.apache.log4j.Logger;

import edu.scripps.yates.nucleome.model.CellCompartment;
import edu.scripps.yates.nucleome.model.CellType;
import edu.scripps.yates.nucleome.model.Experiment;

/**
 * This is the second scoring function proposed by LArry. It is based on the
 * NSAF difference but not the ratio.
 * 
 * @author Salva
 *
 */
public class NSAFScoringFunctionByFoldChangeValue extends ScoringFunction {
	private final static Logger log = Logger.getLogger(NSAFScoringFunctionByFoldChangeValue.class);
	private final _4DNucleomeAnalyzer analyzer;

	public NSAFScoringFunctionByFoldChangeValue(_4DNucleomeAnalyzer analyzer) {
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
			double experimentScore = 0.0;
			int numComparisons = 0;
			boolean ratioInThisExperiment = false;

			for (CellCompartment denominatorCellCompartment : CellCompartment.values()) {
				if (denominatorCellCompartment == Constants.cellCompartmentToStudy) {
					continue;
				}
				double nsafToStudy = experiment.getAvgNSAF(proteinAccessions, Constants.cellCompartmentToStudy, true);
				double nsafDenominator = experiment.getAvgNSAF(proteinAccessions, denominatorCellCompartment, true);
				if (nsafToStudy > 0 || nsafDenominator > 0) {
					ratioInThisExperiment = true;
					numComparisons++;
				} else if (nsafToStudy == 0 && nsafDenominator == 0) {
					continue;
				}
				if (nsafToStudy >= nsafDenominator) {
					experimentScore += 1;
					if (nsafToStudy >= 2 * nsafDenominator && nsafToStudy > 2) {
						experimentScore += 1;
					}
				} else {
					if (Constants.includeNegativeScoring) {
						experimentScore -= 1;
						if (2 * nsafToStudy <= nsafDenominator && nsafDenominator > 2) {
							experimentScore -= 1;
						}
					}
				}
			}
			// normalize by the number of comparisons
			if (numComparisons > 0) {
				experimentScore = experimentScore / numComparisons;
			}
			score += experimentScore;
			if (ratioInThisExperiment) {
				numExperiments++;
			}

		}
		// if we are looking in a particular cellType, divide
		score = score / numExperiments;
		return score;
	}
}
