package edu.scripps.yates.nucleome;

import java.io.File;
import java.util.Collection;
import java.util.regex.Pattern;

import edu.scripps.yates.annotations.uniprot.UniprotProteinRetriever;
import edu.scripps.yates.nucleome.model.CellCompartment;

public class Constants {
	public static int MIN_AVG_SPC = 3;
	public static int MIN_PEPTIDES_PER_PROTEIN = 1;
	public static boolean GO_FILTER;
	public static boolean KERATIN_FILTER;
	public static boolean includeNegativeScoring = false;
	public final static UniprotProteinRetriever upr = new UniprotProteinRetriever(null,
			new File("z:\\share\\Salva\\data\\uniprotKB"), true);
	public final static String decoy = "Reverse|contaminant";
	public static Pattern pattern;
	public static final String CONTROL_FILE = "z:\\share\\Salva\\data\\4D_Nucleome\\NE50 updated 2017_11_21.txt";
	public static final String SEPARATOR = " | ";
	public static boolean writeCoverageFile = false;
	public static boolean writeCombinedDistribution = false;
	public static boolean compareScores = true;
	public static boolean printScoreDistributions = true;
	public static int MAX_TEST_PROTEINS = 200000;
	public static int MIN_TOTAL_SPC = 5;
	public static int MIN_PSM_PER_PROTEIN = 2;
	public static boolean TESTING;

	public static CellCompartment cellCompartmentToStudy = CellCompartment.NE;
	/**
	 * Defines which is the threshold used to consider a protein to be enriched.
	 * <br>
	 * Note that only enriched proteins are considered in the comparisons
	 */
	public static Double ENRICHMENT_SCORE_THRESHOLD;
	public static String DATASET_PATHS_FILE;
	public static String[] geneFilter;

	public static boolean isDecoy(String rawAcc) {
		if (pattern == null) {
			pattern = Pattern.compile(decoy);
		}
		return pattern.matcher(rawAcc).find();
	}

	public static boolean isDecoy(Collection<String> accs) {
		if (pattern == null) {
			pattern = Pattern.compile(decoy);
		}
		for (String acc : accs) {
			boolean isDecoy = pattern.matcher(acc).find();
			if (isDecoy) {
				return true;
			}
		}
		return false;
	}

	public static String formatInfinity(double value) {
		if (!Double.isInfinite(value)) {
			return String.valueOf(value);
		}
		if (Double.POSITIVE_INFINITY == value) {
			return "100";
		} else if (Double.NEGATIVE_INFINITY == value) {
			return "-100";
		}
		return String.valueOf(value);
	}
}
