package edu.scripps.yates.nucleome.filters;

import edu.scripps.yates.utilities.proteomicsmodel.Protein;

public class PSMPerProtein implements Filter {

	private final String name;

	public PSMPerProtein(String runid) {
		name = runid;
	}

	@Override
	public boolean isValid(Protein protein) {
		return protein.getPSMs().size() >= edu.scripps.yates.nucleome.Constants.MIN_PSM_PER_PROTEIN;
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
