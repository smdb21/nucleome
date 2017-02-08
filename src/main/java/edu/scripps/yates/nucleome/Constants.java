package edu.scripps.yates.nucleome;

import java.io.File;
import java.util.Set;
import java.util.regex.Pattern;

import edu.scripps.yates.annotations.uniprot.UniprotProteinRetriever;
import edu.scripps.yates.nucleome.model.CellCompartment;

public class Constants {
	public static int MIN_AVG_SPC = 3;
	public static int MIN_PEPTIDES_PER_PROTEIN = 1;
	public static boolean GO_FILTER;
	public static boolean includeNegativeScoring = false;
	public final static UniprotProteinRetriever upr = new UniprotProteinRetriever(null,
			new File("C:\\Users\\Salva\\Desktop\\tmp\\uniprotKB"), true);
	public final static Pattern decoy = Pattern.compile("Reverse");

	public static int MIN_PSM_PER_PROTEIN = 2;
	public static boolean TESTING;

	public static CellCompartment cellCompartmentToStudy = CellCompartment.NE;
	/**
	 * Defines which is the threshold used to consider a protein to be
	 * enriched.<br>
	 * Note that only enriched proteins are considered in the comparisons
	 */
	public static int ENRICHMENT_SCORE_THRESHOLD;
	public static String DATASET_PATHS_FILE;
	public static boolean USE_GROUPS = false;

	public static boolean isDecoy(String rawAcc) {
		return decoy.matcher(rawAcc).find();
	}

	public static boolean isDecoy(Set<String> accs) {
		for (String acc : accs) {
			boolean isDecoy = decoy.matcher(acc).find();
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
