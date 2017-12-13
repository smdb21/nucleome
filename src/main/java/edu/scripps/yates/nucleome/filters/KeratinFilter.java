package edu.scripps.yates.nucleome.filters;

import java.util.HashSet;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.utilities.proteomicsmodel.Protein;

public class KeratinFilter implements Filter {
	private final static Logger log = Logger.getLogger(KeratinFilter.class);

	private final Set<String> filteredOut = new HashSet<String>();
	private final String name;

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("Kerating filter");
		return sb.toString();
	}

	public KeratinFilter(String name) {
		this.name = name;
	}

	@Override
	public boolean isValid(String accession) {
		if (accession.toLowerCase().contains("keratin")) {
			return false;
		}
		return true;
	}

	@Override
	public boolean isValid(Protein protein) {
		return isValid(protein.getPrimaryAccession().getDescription());
	}

}
