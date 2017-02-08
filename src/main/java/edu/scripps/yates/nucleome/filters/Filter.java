package edu.scripps.yates.nucleome.filters;

import edu.scripps.yates.utilities.proteomicsmodel.Protein;

public interface Filter {
	public boolean isValid(Protein protein);

	public boolean isValid(String acc);
}
