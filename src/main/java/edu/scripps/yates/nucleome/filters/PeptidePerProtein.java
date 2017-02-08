package edu.scripps.yates.nucleome.filters;

import edu.scripps.yates.utilities.proteomicsmodel.Protein;

public class PeptidePerProtein implements Filter {

	private final String name;

	public PeptidePerProtein(String runid) {
		name = runid;
	}

	@Override
	public boolean isValid(Protein protein) {
		return protein.getPeptides().size() >= edu.scripps.yates.nucleome.Constants.MIN_PEPTIDES_PER_PROTEIN;
	}

	/**
	 * @return the name
	 */
	public String getName() {
		return name;
	}

	@Override
	public boolean isValid(String acc) {
		// TODO Auto-generated method stub
		return false;
	}

}
