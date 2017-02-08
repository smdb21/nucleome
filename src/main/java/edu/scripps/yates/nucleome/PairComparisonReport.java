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

import edu.scripps.yates.nucleome.model.CellType;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.proteomicsmodel.Protein;
import edu.scripps.yates.utilities.util.Pair;

public class PairComparisonReport {
	protected final Map<String, Double> scoreMap1;
	protected final Map<String, Double> scoreMap2;
	protected boolean ready = false;
	protected final CellType cellType1;
	protected final CellType cellType2;
	protected List<Pair<String, Double>> uniquesTo1;
	protected List<Pair<String, Double>> uniquesTo2;
	protected List<Pair<String, String>> intersection;

	public PairComparisonReport(CellType cellType1, Map<String, Double> scoreMap1, CellType cellType2,
			Map<String, Double> scoreMap2) {
		this.scoreMap1 = scoreMap1;
		this.scoreMap2 = scoreMap2;
		this.cellType1 = cellType1;
		this.cellType2 = cellType2;
	}

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
			fw.write(uniquesTo1.size() + " proteins enriched in " + Constants.cellCompartmentToStudy + " in "
					+ cellType1 + " and not in " + cellType2 + "\n");
			fw.write(uniquesTo2.size() + " proteins enriched in " + Constants.cellCompartmentToStudy + " in "
					+ cellType2 + " and not in " + cellType1 + "\n");
			fw.write(intersection.size() + " proteins enriched in " + Constants.cellCompartmentToStudy + " in both "
					+ cellType1 + " and " + cellType2 + "\n\n");
			fw.write("\n" + uniquesTo1.size() + " proteins enriched in " + Constants.cellCompartmentToStudy + " in "
					+ cellType1 + " and not in " + cellType2 + ":\n");
			fw.write("#\tACC\tEnrichment_Score\tprotein description\n");
			writeListOfPairs(fw, uniquesTo1);

			fw.write("\n" + uniquesTo2.size() + " proteins enriched in " + Constants.cellCompartmentToStudy + " in "
					+ cellType2 + " and not in " + cellType1 + ":\n");
			fw.write("#\tACC\tEnrichment_Score\tprotein description\n");
			writeListOfPairs(fw, uniquesTo2);

			fw.write("\nProteins enriched in " + Constants.cellCompartmentToStudy + " in both " + cellType1 + " and in "
					+ cellType2 + ":\n");
			fw.write("#\tACC\tDelta_Enrichment_Score (" + cellType1 + "-" + cellType2 + ")\tprotein description\n");
			writeListOfPairs2(fw, intersection);

		} finally {
			fw.close();
		}

	}

	protected void writeListOfPairs(FileWriter fw, List<Pair<String, Double>> listOfPairs) throws IOException {
		int i = 1;
		for (Pair<String, Double> pair : listOfPairs) {
			final String rawAcc = pair.getFirstelement();
			String acc = FastaParser.getACC(rawAcc).getFirstelement();
			final Map<String, Protein> annotatedProtein = Constants.upr.getAnnotatedProtein(acc);
			String proteinName = "";
			if (annotatedProtein.containsKey(acc) && annotatedProtein.get(acc) != null) {
				if (annotatedProtein.get(acc).getPrimaryAccession() != null
						&& annotatedProtein.get(acc).getPrimaryAccession().getDescription() != null) {
					proteinName = annotatedProtein.get(acc).getPrimaryAccession().getDescription();
				}
			}
			fw.write(i++ + "\t" + rawAcc + "\t" + pair.getSecondElement() + "\t" + proteinName + "\n");
		}
	}

	private void compare() {
		uniquesTo1 = uniques(scoreMap1, scoreMap2);
		uniquesTo2 = uniques(scoreMap2, scoreMap1);
		intersection = intersection(scoreMap1, scoreMap2);
		ready = true;
	}

	private List<Pair<String, String>> intersection(Map<String, Double> scoreMap12, Map<String, Double> scoreMap22) {
		List<Pair<String, String>> ret = new ArrayList<Pair<String, String>>();
		for (String proteinAcc : scoreMap12.keySet()) {
			final Double enrichmentScore1 = scoreMap12.get(proteinAcc);
			if (Double.isNaN(enrichmentScore1)) {
				continue;
			}
			// if it is enriched in 1
			if (enrichmentScore1 > Constants.ENRICHMENT_SCORE_THRESHOLD) {
				final Double enrichmentScore2 = scoreMap22.get(proteinAcc);
				if (Double.isNaN(enrichmentScore2)) {
					continue;
				}
				// if it is enriched in 2 also
				if (enrichmentScore2 > Constants.ENRICHMENT_SCORE_THRESHOLD) {
					String scores = String.valueOf(enrichmentScore1) + " | " + String.valueOf(enrichmentScore2);

					Pair<String, String> pair = new Pair<String, String>(proteinAcc, scores);
					ret.add(pair);
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

	private List<Pair<String, Double>> uniques(Map<String, Double> scoreMap1, Map<String, Double> scoreMap2) {
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

	protected double getSumScores(String secondElement) {
		double ret = 0.0;
		final String[] split = secondElement.split("\\|");
		for (String string : split) {
			ret += Double.valueOf(string.trim());
		}
		return ret;
	}

	protected void writeListOfPairs2(FileWriter fw, List<Pair<String, String>> listOfPairs) throws IOException {
		int i = 1;
		for (Pair<String, String> pair : listOfPairs) {
			final String rawAcc = pair.getFirstelement();
			String acc = FastaParser.getACC(rawAcc).getFirstelement();
			final Map<String, Protein> annotatedProtein = Constants.upr.getAnnotatedProtein(acc);
			String proteinName = "";
			if (annotatedProtein.containsKey(acc) && annotatedProtein.get(acc) != null) {
				if (annotatedProtein.get(acc).getPrimaryAccession() != null
						&& annotatedProtein.get(acc).getPrimaryAccession().getDescription() != null) {
					proteinName = annotatedProtein.get(acc).getPrimaryAccession().getDescription();
				}
			}
			fw.write(i++ + "\t" + rawAcc + "\t" + pair.getSecondElement() + "\t" + proteinName + "\n");
		}
	}
}
