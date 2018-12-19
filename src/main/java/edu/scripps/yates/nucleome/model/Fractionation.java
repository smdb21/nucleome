package edu.scripps.yates.nucleome.model;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.dtaselectparser.DTASelectParser;
import edu.scripps.yates.nucleome.Constants;
import edu.scripps.yates.nucleome.filters.Filter;
import edu.scripps.yates.nucleome.filters.PSMPerProtein;
import edu.scripps.yates.nucleome.filters.PeptidePerProtein;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.proteomicsmodel.Amount;
import edu.scripps.yates.utilities.proteomicsmodel.PSM;
import edu.scripps.yates.utilities.proteomicsmodel.Peptide;
import edu.scripps.yates.utilities.proteomicsmodel.Protein;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AmountType;
import edu.scripps.yates.utilities.remote.RemoteSSHFileReference;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.set.hash.THashSet;

public class Fractionation {
	private final CellCompartment cellCompartment;
	private final DTASelectParser parser;
	private HashSet<String> proteinAccs;
	private final CellType cellType;

	private final List<Filter> filters;
	private final String name;
	private Map<String, String> proteinSequences;
	private final Wash wash;
	private List<Protein> proteins;

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
		parser = new DTASelectParser(name, remote);
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

	private Map<String, String> getProteinSequences() throws IOException {
		if (proteinSequences == null) {
			proteinSequences = Constants.upr.getAnnotatedProteinSequence(parser.getProteinMap().keySet());
		}
		return proteinSequences;
	}

	private List<Filter> getFilters() {

		final List<Filter> filters = new ArrayList<Filter>();
		final PeptidePerProtein peptidesPerProteinFilter = new PeptidePerProtein("ASDF2");
		filters.add(peptidesPerProteinFilter);
		final PSMPerProtein psmsPerProteinFilter = new PSMPerProtein("ASDF3");
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
	public DTASelectParser getParser() {
		return parser;
	}

	public Set<String> getProteinAccs() throws IOException {
		if (proteinAccs == null) {
			proteinAccs = new HashSet<String>();
			final Set<String> keySet = parser.getProteinMap().keySet();

			proteinAccs.addAll(keySet);

		}
		return proteinAccs;
	}

	public List<Protein> getProteins() throws IOException {
		if (proteins == null) {
			proteins = new ArrayList<Protein>();

			for (final Protein protein : parser.getProteins()) {
				if (protein.getLength() == null || protein.getLength() <= 0) {
					if (getProteinSequences().containsKey(protein.getAccession())) {
						protein.setLength(getProteinSequences().get(protein.getAccession()).length());
					}
				}
				proteins.add(protein);

			}
		}
		return proteins;
	}

	public Set<PSM> getPSMs(String proteinAcc, boolean skipFilters) throws IOException {
		final Set<PSM> psms = new HashSet<PSM>();

		if (getProteinAccs().contains(proteinAcc)) {
			final Protein protein = parser.getProteinMap().get(proteinAcc);
			if (!skipFilters) {
				for (final Filter filter : filters) {
					if (!filter.isValid(protein)) {
						return Collections.emptySet();
					}
				}

				psms.addAll(protein.getPSMs());
			}

		}
		return psms;
	}

	public int getPeptideCount(String proteinAcc, boolean skipFilters) throws IOException {
		final Set<String> peptideSequences = new HashSet<String>();

		if (getProteinAccs().contains(proteinAcc)) {
			final Protein protein = parser.getProteinMap().get(proteinAcc);
			if (!skipFilters) {
				for (final Filter filter : filters) {
					if (!filter.isValid(protein)) {
						return 0;
					}
				}
			}
			for (final Peptide peptide : protein.getPeptides()) {
				peptideSequences.add(peptide.getSequence());
			}

		}
		return peptideSequences.size();
	}

	public int getSpectralCount(Collection<String> proteinAccessions, boolean skipFilters) throws IOException {
		final Set<PSM> psms = new HashSet<PSM>();
		final Set<String> accs = new HashSet<String>();
		for (final String proteinAccession : proteinAccessions) {
			if (accs.contains(proteinAccession)) {
				continue;
			}
			if (getProteinAccs().contains(proteinAccession)) {
				accs.add(proteinAccession);
				final Protein protein = parser.getProteinMap().get(proteinAccession);
				boolean valid = true;
				if (!skipFilters) {
					for (final Filter filter : filters) {
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
		return psms.size();
	}

	public double getSumNSAF(Collection<String> proteinAccessions, boolean skipFilters) throws IOException {
		double ret = 0.0;
		final Set<String> accs = new THashSet<String>();
		for (final String proteinAccession : proteinAccessions) {
			if (accs.contains(proteinAccession)) {
				continue;
			}
			if (getProteinAccs().contains(proteinAccession)) {
				accs.add(proteinAccession);
				final Protein protein = parser.getProteinMap().get(proteinAccession);
				// we make the average of them because DTASelect gives us a
				// protein instance per MSRUN and all of them share the same
				// NSAF. To sum all up would be redundant.
				final TDoubleArrayList toAverage = new TDoubleArrayList();
				boolean valid = true;
				if (!skipFilters) {
					for (final Filter filter : filters) {
						if (!filter.isValid(protein)) {
							valid = false;
						}
					}
				}
				if (valid) {
					final Double nsafFromProtein = getNSAFFromProtein(protein);
					if (nsafFromProtein != null) {
						toAverage.add(nsafFromProtein);
					}
				}

				ret += Maths.mean(toAverage);
			}
		}
		return ret;
	}

	public double getAverageNSAF(Collection<String> proteinAccessions, boolean skipFilters) throws IOException {
		final List<Double> nsafs = new ArrayList<Double>();
		final Set<String> accs = new HashSet<String>();
		for (final String proteinAccession : proteinAccessions) {
			if (accs.contains(proteinAccession)) {
				continue;
			}
			if (getProteinAccs().contains(proteinAccession)) {
				accs.add(proteinAccession);
				final Protein protein = parser.getProteinMap().get(proteinAccession);
				// System.out.println(protein);
				boolean valid = true;
				if (!skipFilters) {
					for (final Filter filter : filters) {
						if (!filter.isValid(protein)) {
							valid = false;
						}
					}
				}
				if (valid) {
					final Double nsafFromProtein = getNSAFFromProtein(protein);
					// System.out.println(nsafFromProtein + "\t" +
					// protein.hashCode());
					if (nsafFromProtein != null) {
						nsafs.add(nsafFromProtein);
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
			for (final Amount amount : protein.getAmounts()) {
				if (amount.getAmountType() == AmountType.NSAF) {
					return amount.getValue();
				}
			}
		}
		return protein.getNsaf().doubleValue();
	}

	public int getPeptideCount(Collection<String> proteinAccessions, boolean skipFilters) throws IOException {
		final Set<String> peptides = new HashSet<String>();
		final Set<String> accs = new HashSet<String>();

		for (final String proteinAccession : proteinAccessions) {
			if (accs.contains(proteinAccession)) {
				continue;
			}
			if (getProteinAccs().contains(proteinAccession)) {
				accs.add(proteinAccession);
				final Protein protein = parser.getProteinMap().get(proteinAccession);
				boolean valid = true;
				if (!skipFilters) {
					for (final Filter filter : filters) {
						if (!filter.isValid(protein)) {
							valid = false;
						}
					}
				}
				if (valid) {
					for (final Peptide peptide : protein.getPeptides()) {
						peptides.add(peptide.getSequence());
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
		final StringBuilder sb = new StringBuilder();
		final Set<String> accs = new HashSet<String>();
		final DecimalFormat formatter = new DecimalFormat("#.#");
		for (final String proteinAccession : proteinAccessions) {
			if (accs.contains(proteinAccession)) {
				continue;
			}
			if (getProteinAccs().contains(proteinAccession)) {
				accs.add(proteinAccession);
				final Protein protein = parser.getProteinMap().get(proteinAccession);
				boolean valid = true;
				if (!skipFilters) {
					for (final Filter filter : filters) {
						if (!filter.isValid(protein)) {
							valid = false;
						}
					}
				}
				if (valid) {

					final Float coverage = protein.getCoverage();
					if (coverage != null) {
						if (!"".equals(sb.toString())) {
							sb.append(",");
						}
						sb.append(formatter.format(coverage) + "%");
						break;
					}
				}

			}
		}

		return sb.toString();

	}
}
