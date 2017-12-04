package edu.scripps.yates.nucleome.model;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import edu.scripps.yates.dtaselect.ProteinDTASelectParser;
import edu.scripps.yates.nucleome.Constants;
import edu.scripps.yates.nucleome.filters.Filter;
import edu.scripps.yates.nucleome.filters.PSMPerProtein;
import edu.scripps.yates.nucleome.filters.PeptidePerProtein;
import edu.scripps.yates.utilities.grouping.GroupableProtein;
import edu.scripps.yates.utilities.grouping.ProteinGroup;
import edu.scripps.yates.utilities.proteomicsmodel.PSM;
import edu.scripps.yates.utilities.proteomicsmodel.Peptide;
import edu.scripps.yates.utilities.proteomicsmodel.Protein;
import edu.scripps.yates.utilities.remote.RemoteSSHFileReference;

public class Fractionation {
	private final CellCompartment cellCompartment;
	private final ProteinDTASelectParser parser;
	private HashSet<String> proteinAccs;
	private final CellType cellType;

	private final List<Filter> filters;
	private final String name;

	public Fractionation(String experimentName, int num, CellCompartment cellCompartment, RemoteSSHFileReference remote,
			CellType cellType) throws IOException {
		this.cellCompartment = cellCompartment;
		name = experimentName + cellCompartment.name() + "_rep" + num;
		parser = new ProteinDTASelectParser(name, remote);
		// removed since we want to see how is the enrichment score of the
		// decoys:
		if (Constants.decoy != null) {
			parser.setDecoyPattern(Constants.decoy);
		}

		this.cellType = cellType;
		filters = getFilters();

	}

	private List<Filter> getFilters() {

		List<Filter> filters = new ArrayList<Filter>();
		PeptidePerProtein peptidesPerProteinFilter = new PeptidePerProtein("ASDF2");
		filters.add(peptidesPerProteinFilter);
		PSMPerProtein psmsPerProteinFilter = new PSMPerProtein("ASDF3");
		filters.add(psmsPerProteinFilter);

		return filters;
	}

	public Fractionation(String experimentName, int num, CellCompartment cellCompartment, File file, CellType cellType)
			throws IOException {
		this(experimentName, num, cellCompartment, new RemoteSSHFileReference(file), cellType);
	}

	/**
	 * @return the cellCompartment
	 */
	public CellCompartment getCellCompartment() {
		return cellCompartment;
	}

	/**
	 * @return the remote
	 */
	public ProteinDTASelectParser getParser() {
		return parser;
	}

	public Set<String> getProteinAccs() throws IOException {
		if (proteinAccs == null) {
			proteinAccs = new HashSet<String>();
			final Set<String> keySet = parser.getProteins().keySet();
			for (String acc : keySet) {
				proteinAccs.add(acc);
			}
		}
		return proteinAccs;
	}

	public Set<Protein> getProteins() throws IOException {
		Set<Protein> ret = new HashSet<Protein>();

		for (Set<Protein> proteinSet : parser.getProteins().values()) {
			for (Protein protein : proteinSet) {
				ret.add(protein);
			}
		}
		return ret;
	}

	public int getSpectralCount(String proteinAcc, boolean skipFilters) throws IOException {
		Set<PSM> psms = new HashSet<PSM>();
		if (proteinAccs.contains(proteinAcc)) {
			final Set<Protein> proteins = parser.getProteins().get(proteinAcc);
			for (Protein protein : proteins) {
				if (!skipFilters) {
					for (Filter filter : filters) {
						if (!filter.isValid(protein)) {
							return 0;
						}
					}
				}
				psms.addAll(protein.getPSMs());
			}

		}
		return psms.size();
	}

	public int getPeptideCount(String proteinAcc, boolean skipFilters) throws IOException {
		Set<String> peptideSequences = new HashSet<String>();
		if (proteinAccs.contains(proteinAcc)) {
			final Set<Protein> proteins = parser.getProteins().get(proteinAcc);
			for (Protein protein : proteins) {
				if (!skipFilters) {
					for (Filter filter : filters) {
						if (!filter.isValid(protein)) {
							return 0;
						}
					}
				}
				for (Peptide peptide : protein.getPeptides()) {
					peptideSequences.add(peptide.getSequence());
				}
			}

		}
		return peptideSequences.size();
	}

	public int getSpectralCount(ProteinGroup proteinGroup, boolean skipFilters) throws IOException {
		Set<PSM> psms = new HashSet<PSM>();
		for (GroupableProtein groupableProtein : proteinGroup) {
			if (proteinAccs.contains(groupableProtein.getAccession())) {
				final Set<Protein> proteins = parser.getProteins().get(groupableProtein.getAccession());
				for (Protein protein : proteins) {
					boolean valid = true;
					if (!skipFilters) {
						for (Filter filter : filters) {
							if (!filter.isValid(protein)) {
								valid = false;
							}
						}
					}
					if (valid) {
						psms.addAll(protein.getPSMs());
					}
				}
			}
		}
		return psms.size();
	}

	public int getPeptideCount(ProteinGroup proteinGroup, boolean skipFilters) throws IOException {
		Set<String> peptides = new HashSet<String>();
		for (GroupableProtein groupableProtein : proteinGroup) {
			if (proteinAccs.contains(groupableProtein.getAccession())) {
				final Set<Protein> proteins = parser.getProteins().get(groupableProtein.getAccession());
				for (Protein protein : proteins) {
					boolean valid = true;
					if (!skipFilters) {
						for (Filter filter : filters) {
							if (!filter.isValid(protein)) {
								valid = false;
							}
						}
					}
					if (valid) {
						for (Peptide peptide : protein.getPeptides()) {
							peptides.add(peptide.getSequence());
						}
					}
				}
			}
		}
		return peptides.size();
	}

	/**
	 * @return the cellType
	 */
	public CellType getCellType() {
		return cellType;
	}

	public String getName() {
		return name;
	}
}
