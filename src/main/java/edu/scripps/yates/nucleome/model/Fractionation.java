package edu.scripps.yates.nucleome.model;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.dtaselect.ProteinDTASelectParser;
import edu.scripps.yates.dtaselect.ProteinImplFromDTASelect;
import edu.scripps.yates.nucleome.Constants;
import edu.scripps.yates.nucleome.filters.Filter;
import edu.scripps.yates.nucleome.filters.PSMPerProtein;
import edu.scripps.yates.nucleome.filters.PeptidePerProtein;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.model.enums.AmountType;
import edu.scripps.yates.utilities.proteomicsmodel.Amount;
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
	private Map<String, String> proteinSequences;
	private final Wash wash;

	public Fractionation(String experimentName, int num, CellCompartment cellCompartment, RemoteSSHFileReference remote,
			CellType cellType) throws IOException {
		this(experimentName, num, cellCompartment, remote, cellType, null);
	}

	public Fractionation(String experimentName, int num, CellCompartment cellCompartment, RemoteSSHFileReference remote,
			CellType cellType, Wash wash) throws IOException {
		this.cellCompartment = cellCompartment;

		if (wash != null) {
			name = experimentName + "_" + cellCompartment.name() + "_rep" + num + "_" + wash;
		} else {
			name = experimentName + "_" + cellCompartment.name() + "_rep" + num;
		}
		parser = new ProteinDTASelectParser(name, remote);
		// parser.setSeparateByFractionationStep(false);
		// removed since we want to see how is the enrichment score of the
		// decoys:
		if (Constants.decoy != null) {
			parser.setDecoyPattern(Constants.decoy);
		}

		this.cellType = cellType;
		this.wash = wash;
		filters = getFilters();

	}

	private Map<String, String> getProteinSequences() {
		if (proteinSequences == null) {
			proteinSequences = Constants.upr.getAnnotatedProteinSequence(parser.getProteins().keySet());
		}
		return this.proteinSequences;
	}

	private List<Filter> getFilters() {

		List<Filter> filters = new ArrayList<Filter>();
		PeptidePerProtein peptidesPerProteinFilter = new PeptidePerProtein("ASDF2");
		filters.add(peptidesPerProteinFilter);
		PSMPerProtein psmsPerProteinFilter = new PSMPerProtein("ASDF3");
		filters.add(psmsPerProteinFilter);

		return filters;
	}

	public Fractionation(String experimentName, int num, CellCompartment cellCompartment, File file, CellType cellType,
			Wash wash) throws IOException {
		this(experimentName, num, cellCompartment, new RemoteSSHFileReference(file), cellType, wash);
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
				if (protein.getLength() <= 0) {
					if (getProteinSequences().containsKey(protein.getAccession())) {
						protein.setLength(getProteinSequences().get(protein.getAccession()).length());
					}
				}
				ret.add(protein);
			}
		}
		return ret;
	}

	public int getSpectralCount(String proteinAcc, boolean skipFilters) throws IOException {
		Set<PSM> psms = new HashSet<PSM>();

		if (getProteinAccs().contains(proteinAcc)) {
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

		if (getProteinAccs().contains(proteinAcc)) {
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

	public int getSpectralCount(Collection<String> proteinAccessions, boolean skipFilters) throws IOException {
		Set<PSM> psms = new HashSet<PSM>();
		Set<String> accs = new HashSet<String>();
		for (String proteinAccession : proteinAccessions) {
			if (accs.contains(proteinAccession)) {
				continue;
			}
			if (getProteinAccs().contains(proteinAccession)) {
				accs.add(proteinAccession);
				final Set<Protein> proteins = parser.getProteins().get(proteinAccession);
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

	public double getSumNSAF(Collection<String> proteinAccessions, boolean skipFilters) throws IOException {
		double ret = 0.0;
		Set<String> accs = new HashSet<String>();
		for (String proteinAccession : proteinAccessions) {
			if (accs.contains(proteinAccession)) {
				continue;
			}
			if (getProteinAccs().contains(proteinAccession)) {
				accs.add(proteinAccession);
				final Set<Protein> proteins = parser.getProteins().get(proteinAccession);
				// we make the average of them because DTASelect gives us a
				// protein instance per MSRUN and all of them share the same
				// NSAF. To sum all up would be redundant.
				List<Double> toAverage = new ArrayList<Double>();
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
						Double nsafFromProtein = getNSAFFromProtein(protein);
						if (nsafFromProtein != null) {
							toAverage.add(nsafFromProtein);
						}
					}
				}
				ret += Maths.mean(toAverage.toArray(new Double[0]));
			}
		}
		return ret;
	}

	public double getAverageNSAF(Collection<String> proteinAccessions, boolean skipFilters) throws IOException {
		List<Double> nsafs = new ArrayList<Double>();
		Set<String> accs = new HashSet<String>();
		for (String proteinAccession : proteinAccessions) {
			if (accs.contains(proteinAccession)) {
				continue;
			}
			if (getProteinAccs().contains(proteinAccession)) {
				accs.add(proteinAccession);
				final Set<Protein> proteins = parser.getProteins().get(proteinAccession);
				for (Protein protein : proteins) {
					// System.out.println(protein);
					boolean valid = true;
					if (!skipFilters) {
						for (Filter filter : filters) {
							if (!filter.isValid(protein)) {
								valid = false;
							}
						}
					}
					if (valid) {
						Double nsafFromProtein = getNSAFFromProtein(protein);
						// System.out.println(nsafFromProtein + "\t" +
						// protein.hashCode());
						if (nsafFromProtein != null) {
							nsafs.add(nsafFromProtein);
						}
					}
				}
			}
		}
		if (!nsafs.isEmpty()) {
			return Maths.mean(nsafs.toArray(new Double[0]));
		}
		return 0.0;
	}

	private Double getNSAFFromProtein(Protein protein) {
		if (protein != null && protein.getAmounts() != null) {
			for (Amount amount : protein.getAmounts()) {
				if (amount.getAmountType() == AmountType.NSAF) {
					return amount.getValue();
				}
			}
		}
		return null;
	}

	public int getPeptideCount(Collection<String> proteinAccessions, boolean skipFilters) throws IOException {
		Set<String> peptides = new HashSet<String>();
		Set<String> accs = new HashSet<String>();

		for (String proteinAccession : proteinAccessions) {
			if (accs.contains(proteinAccession)) {
				continue;
			}
			if (getProteinAccs().contains(proteinAccession)) {
				accs.add(proteinAccession);
				final Set<Protein> proteins = parser.getProteins().get(proteinAccession);
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

	public Wash getWash() {
		return wash;
	}

	public String getCoverage(Collection<String> proteinAccessions, boolean skipFilters) throws IOException {
		StringBuilder sb = new StringBuilder();
		Set<String> accs = new HashSet<String>();
		DecimalFormat formatter = new DecimalFormat("#.#");
		for (String proteinAccession : proteinAccessions) {
			if (accs.contains(proteinAccession)) {
				continue;
			}
			if (getProteinAccs().contains(proteinAccession)) {
				accs.add(proteinAccession);
				final Set<Protein> proteins = parser.getProteins().get(proteinAccession);
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
						if (protein instanceof ProteinImplFromDTASelect) {
							ProteinImplFromDTASelect proteinDTASelect = (ProteinImplFromDTASelect) protein;
							double coverage = proteinDTASelect.getCoverage();
							if (!"".equals(sb.toString())) {
								sb.append(",");
							}
							sb.append(formatter.format(coverage) + "%");
							break;
						}
					}
				}
			}
		}

		return sb.toString();

	}
}
