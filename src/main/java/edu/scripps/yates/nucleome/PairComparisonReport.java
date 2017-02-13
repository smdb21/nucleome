package edu.scripps.yates.nucleome;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.nucleome.model.CellType;
import edu.scripps.yates.utilities.util.Pair;
import edu.scripps.yates.utilities.venndata.VennData;

public class PairComparisonReport {
	private final _4DNucleomeAnalyzer nucleomeAnalyzer;
	protected final Map<String, Double> scoreMap1;
	protected final Map<String, Double> scoreMap2;
	protected boolean ready = false;
	protected final CellType cellType1;
	protected final CellType cellType2;
	protected List<Pair<String, Double>> uniquesTo1;
	protected List<Pair<String, Double>> uniquesTo2;
	protected List<Pair<String, String>> intersection;
	protected VennData venn;

	public PairComparisonReport(_4DNucleomeAnalyzer nucleomeAnalyzer, CellType cellType1, Map<String, Double> scoreMap1,
			CellType cellType2, Map<String, Double> scoreMap2) {
		this.nucleomeAnalyzer = nucleomeAnalyzer;
		this.scoreMap1 = scoreMap1;
		this.scoreMap2 = scoreMap2;
		this.cellType1 = cellType1;
		this.cellType2 = cellType2;
	}

	public VennData getVennData() throws IOException {
		if (!ready) {
			compare();
		}
		if (venn == null) {
			Set<String> enriched1 = getEnriched(scoreMap1);
			Set<String> enriched2 = getEnriched(scoreMap2);
			venn = new VennData(cellType1.name() + " vs " + cellType2.name(), cellType1.name(), enriched1,
					cellType2.name(), enriched2, null, null);
		}
		return venn;
	}

	protected Set<String> getEnriched(Map<String, Double> scoreMap) throws IOException {
		Set<String> ret = new HashSet<String>();
		for (String proteinAcc : scoreMap.keySet()) {
			Double enrichmentScore = scoreMap.get(proteinAcc);
			if (!this.nucleomeAnalyzer.isValid(proteinAcc)) {
				continue;
			}
			if (!Double.isNaN(enrichmentScore) && enrichmentScore >= Constants.ENRICHMENT_SCORE_THRESHOLD) {
				ret.add(proteinAcc);
			}
		}
		return ret;
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
			fw.write(getEnriched(scoreMap1).size() + " proteins enriched in " + Constants.cellCompartmentToStudy
					+ " in " + cellType1 + "\n");
			fw.write(getEnriched(scoreMap2).size() + " proteins enriched in " + Constants.cellCompartmentToStudy
					+ " in " + cellType2 + "\n");
			fw.write(uniquesTo1.size() + " proteins enriched in " + Constants.cellCompartmentToStudy + " in "
					+ cellType1 + " and not in " + cellType2 + "\n");
			fw.write(uniquesTo2.size() + " proteins enriched in " + Constants.cellCompartmentToStudy + " in "
					+ cellType2 + " and not in " + cellType1 + "\n");
			fw.write(intersection.size() + " proteins enriched in " + Constants.cellCompartmentToStudy + " in both "
					+ cellType1 + " and " + cellType2 + "\n\n");
			fw.write("\n" + uniquesTo1.size() + " proteins enriched in " + Constants.cellCompartmentToStudy + " in "
					+ cellType1 + " and not in " + cellType2 + ":\n");

			fw.write("\nCopy and paste this URL to see the graphical representation of the comparison:\n"
					+ getVennData().getImageURL().toString() + "\n\n");

			fw.write("#\tACC(s)\tEnrichment_Score\tgene name(s)\tprotein description(s)\n");
			writeListOfPairs(fw, uniquesTo1);

			fw.write("\n" + uniquesTo2.size() + " proteins enriched in " + Constants.cellCompartmentToStudy + " in "
					+ cellType2 + " and not in " + cellType1 + ":\n");
			fw.write("#\tACC\tEnrichment_Score\tprotein description\n");
			writeListOfPairs(fw, uniquesTo2);

			fw.write("\n" + intersection.size() + " proteins enriched in " + Constants.cellCompartmentToStudy
					+ " in both " + cellType1 + " and in " + cellType2 + ":\n");
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
			String proteinNameString = nucleomeAnalyzer.getProteinNameString(rawAcc);
			String geneNameString = nucleomeAnalyzer.getGeneNameString(rawAcc);
			fw.write(i++ + "\t" + rawAcc + "\t" + pair.getSecondElement() + "\t" + geneNameString + "\t"
					+ proteinNameString + "\n");
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
			if (enrichmentScore1 >= Constants.ENRICHMENT_SCORE_THRESHOLD) {
				final Double enrichmentScore2 = scoreMap22.get(proteinAcc);
				if (Double.isNaN(enrichmentScore2)) {
					continue;
				}
				// if it is enriched in 2 also
				if (enrichmentScore2 >= Constants.ENRICHMENT_SCORE_THRESHOLD) {
					String scores = String.valueOf(enrichmentScore1) + Constants.SEPARATOR
							+ String.valueOf(enrichmentScore2);

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
			if (score1 >= Constants.ENRICHMENT_SCORE_THRESHOLD) {
				final Double score2 = scoreMap2.get(proteinAcc);
				// if is enriched in 2, discarded
				if (!Double.isNaN(score2) && score2 >= Constants.ENRICHMENT_SCORE_THRESHOLD) {
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
			String proteinNameString = nucleomeAnalyzer.getProteinNameString(rawAcc);
			String geneNameString = nucleomeAnalyzer.getGeneNameString(rawAcc);
			fw.write(i++ + "\t" + rawAcc + "\t" + pair.getSecondElement() + "\t" + geneNameString + "\t"
					+ proteinNameString + "\n");
		}
	}
}
