package edu.scripps.yates.nucleome;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.nucleome.model.CellType;
import edu.scripps.yates.utilities.util.Pair;
import edu.scripps.yates.utilities.venndata.VennData;

public class TripleComparisonReport extends PairComparisonReport {

	private List<Pair<String, Double>> uniquesTo3;
	private final Map<String, Double> scoreMap3;
	private final CellType cellType3;

	public TripleComparisonReport(_4DNucleomeAnalyzer nucleomeAnalyzer, CellType cellType1,
			Map<String, Double> scoreMap1, CellType cellType2, Map<String, Double> scoreMap2, CellType cellType3,
			Map<String, Double> scoreMap3) {
		super(nucleomeAnalyzer, cellType1, scoreMap1, cellType2, scoreMap2);

		this.scoreMap3 = scoreMap3;
		this.cellType3 = cellType3;
	}

	public VennData getVennData() throws IOException {
		if (!ready) {
			compare();
		}
		if (venn == null) {
			Set<String> enriched1 = getEnriched(scoreMap1);
			Set<String> enriched2 = getEnriched(scoreMap2);
			Set<String> enriched3 = getEnriched(scoreMap3);
			venn = new VennData(cellType1.name() + " vs " + cellType2.name() + " vs " + cellType3.name(),
					cellType1.name(), enriched1, cellType2.name(), enriched2, cellType3.name(), enriched3);
		}
		return venn;
	}

	@Override
	public void printToFile(File file) throws IOException {
		if (!ready) {
			compare();
		}
		// pritn to the file
		FileWriter fw = null;
		try {
			fw = new FileWriter(file);
			fw.write("Report made on " + new Date() + "\n");
			fw.write("We consider enriched proteins when having an enrichment score > "
					+ Constants.ENRICHMENT_SCORE_THRESHOLD + "\n\n");
			fw.write(getEnriched(scoreMap1).size() + " proteins enriched in " + Constants.cellCompartmentToStudy
					+ " in " + cellType1 + "\n");
			fw.write(getEnriched(scoreMap2).size() + " proteins enriched in " + Constants.cellCompartmentToStudy
					+ " in " + cellType2 + "\n");
			fw.write(getEnriched(scoreMap3).size() + " proteins enriched in " + Constants.cellCompartmentToStudy
					+ " in " + cellType3 + "\n");

			fw.write(uniquesTo1.size() + " proteins enriched in " + Constants.cellCompartmentToStudy + " in "
					+ cellType1 + " and not in " + cellType2 + " and not in " + cellType3 + "\n");
			fw.write(uniquesTo2.size() + " proteins enriched in " + Constants.cellCompartmentToStudy + " in "
					+ cellType2 + " and not in " + cellType1 + " and not in " + cellType3 + "\n");
			fw.write(uniquesTo3.size() + " proteins enriched in " + Constants.cellCompartmentToStudy + " in "
					+ cellType3 + " and not in " + cellType1 + " and not in " + cellType2 + "\n");
			fw.write(intersection.size() + " proteins enriched in " + Constants.cellCompartmentToStudy + " in all "
					+ cellType1 + " and " + cellType2 + " and " + cellType3 + "\n\n");

			fw.write("\nCopy and paste this URL to see the graphical representation of the comparison:\n"
					+ getVennData().getImageURL().toString() + "\n\n");
			fw.write("#\tACC\tEnrichment_Score\tgene name(s)\tprotein description(s)\n");
			writeListOfPairs(fw, uniquesTo1);

			fw.write("\n" + uniquesTo2.size() + " proteins enriched in " + Constants.cellCompartmentToStudy + " in "
					+ cellType2 + " and not in " + cellType1 + " and not in " + cellType3 + ":\n");
			fw.write("#\tACC\tEnrichment_Score\tprotein description\n");
			writeListOfPairs(fw, uniquesTo2);

			fw.write("\n" + uniquesTo3.size() + " proteins enriched in " + Constants.cellCompartmentToStudy + " in "
					+ cellType3 + " and not in " + cellType1 + " and not in " + cellType2 + ":\n");
			fw.write("#\tACC\tEnrichment_Score\tprotein description\n");
			writeListOfPairs(fw, uniquesTo3);

			fw.write("\nProteins enriched in " + Constants.cellCompartmentToStudy + " in all " + cellType1 + " and in "
					+ cellType2 + " and in " + cellType3 + ":\n");
			fw.write("#\tACC\tDelta_Enrichment_Score (" + cellType1 + "-" + cellType2 + ")\tprotein description\n");
			writeListOfPairs2(fw, intersection);

		} finally {
			fw.close();
		}

	}

	private void compare() {
		uniquesTo1 = uniques(scoreMap1, scoreMap2, scoreMap3);
		uniquesTo2 = uniques(scoreMap2, scoreMap1, scoreMap3);
		uniquesTo3 = uniques(scoreMap3, scoreMap1, scoreMap2);
		intersection = intersection(scoreMap1, scoreMap2, scoreMap3);
		ready = true;
	}

	private List<Pair<String, String>> intersection(Map<String, Double> scoreMap1, Map<String, Double> scoreMap2,
			Map<String, Double> scoreMap3) {
		List<Pair<String, String>> ret = new ArrayList<Pair<String, String>>();
		for (String proteinAcc : scoreMap1.keySet()) {
			final Double enrichmentScore1 = scoreMap1.get(proteinAcc);
			if (Double.isNaN(enrichmentScore1)) {
				continue;
			}
			// if it is enriched in 1
			if (enrichmentScore1 > Constants.ENRICHMENT_SCORE_THRESHOLD) {
				final Double enrichmentScore2 = scoreMap2.get(proteinAcc);
				if (Double.isNaN(enrichmentScore2)) {
					continue;
				}
				// if it is enriched in 2 also
				if (enrichmentScore2 > Constants.ENRICHMENT_SCORE_THRESHOLD) {
					final Double enrichmentScore3 = scoreMap3.get(proteinAcc);
					if (Double.isNaN(enrichmentScore3)) {
						continue;
					}
					// if it is enriched in 3 also
					if (enrichmentScore3 > Constants.ENRICHMENT_SCORE_THRESHOLD) {
						String scores = String.valueOf(enrichmentScore1) + " | " + String.valueOf(enrichmentScore2)
								+ " | " + String.valueOf(enrichmentScore3);
						Pair<String, String> pair = new Pair<String, String>(proteinAcc, scores);
						ret.add(pair);
					}

				}
			}
		}
		// sort by the delta, from higher to lower
		Collections.sort(ret, new Comparator<Pair<String, String>>() {

			@Override
			public int compare(Pair<String, String> o1, Pair<String, String> o2) {
				return Double.compare(getSumScores(o2.getSecondElement()), getSumScores(o1.getSecondElement()));
			}
		});
		return ret;
	}

	private List<Pair<String, Double>> uniques(Map<String, Double> scoreMap1, Map<String, Double> scoreMap2,
			Map<String, Double> scoreMap3) {
		List<Pair<String, Double>> ret = new ArrayList<Pair<String, Double>>();
		for (String proteinAcc : scoreMap1.keySet()) {
			final Double score1 = scoreMap1.get(proteinAcc);
			if (Double.isNaN(score1)) {
				continue;
			}
			// if it is enriched in 1
			if (score1 > Constants.ENRICHMENT_SCORE_THRESHOLD) {
				final Double score2 = scoreMap2.get(proteinAcc);
				// if is enriched in 2, discarded
				if (!Double.isNaN(score2) && score2 > Constants.ENRICHMENT_SCORE_THRESHOLD) {
					continue;
				}
				final Double score3 = scoreMap3.get(proteinAcc);
				// if is enriched in 3, discarded
				if (!Double.isNaN(score3) && score3 > Constants.ENRICHMENT_SCORE_THRESHOLD) {
					continue;
				}
				Pair<String, Double> pair = new Pair<String, Double>(proteinAcc, score1);
				ret.add(pair);
			}
		}
		// sort by the score
		Collections.sort(ret, new Comparator<Pair<String, Double>>() {

			@Override
			public int compare(Pair<String, Double> o1, Pair<String, Double> o2) {
				return Double.compare(o2.getSecondElement(), o1.getSecondElement());
			}
		});
		return ret;
	}

}
