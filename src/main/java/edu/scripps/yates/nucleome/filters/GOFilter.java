package edu.scripps.yates.nucleome.filters;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.go.GORetriever;
import edu.scripps.yates.annotations.go.GoEntry;
import edu.scripps.yates.nucleome.Constants;
import edu.scripps.yates.utilities.annotations.uniprot.xml.DbReferenceType;
import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.proteomicsmodel.Protein;
import psidev.psi.tools.ontology_manager.interfaces.OntologyTermI;

public class GOFilter implements Filter {
	private final static Logger log = Logger.getLogger(GOFilter.class);
	private final static String[] GOToInclude = { "GO:0005635", "GO:0031965", "GO:0090286", "GO:0006998", "GO:0090292",
			"GO:0005639", "GO:0005638", "GO:0016363", "GO:0048471", "GO:0006997", "GO:1900180", "GO:0035105",
			"GO:0044614", "GO:0005637", "GO:0005643", "GO:0017056", "GO:0097726" };

	private final static String[] GOToExclude = {
			// "GO:0035097",
			// "GO:0034708",
			// "GO:1990234",
			"GO:0016591", "GO:0055029", "GO:0000428", "GO:0030880", "GO:0031248", "GO:1902493",
			// "GO:0016605",
			"GO:0005667",
			// "GO:0005730",
			"GO:0005615", "GO:0044380", "GO:0031012", "GO:0000176", "GO:0101019", "GO:0070062", "GO:1990563",
			"GO:0000178", "GO:0000177",
			// "GO:0005886",
			"GO:0010008", "GO:0005681",
			// "GO:0016604",
			"GO:0005814",
			// "GO:0006397",
			"GO:0005576", "GO:0035327", "GO:0006397",
			// "GO:0008033",
			"GO:0008380", "GO:0042995", "GO:0005911", "GO:0006487", "GO:0006508", "GO:0006412",
			// "GO:0042254",
			"GO:0005925", "GO:0016192", "GO:0006281", "GO:0044822", "GO:0014704", "GO:0006355", "GO:0070062",
			"GO:0005576",
			// new ones added on 6 Mar 2017
			"GO:0005746", // mitochondrial respiratory chain
			"GO:0005743", // mitochondrial inner membrane
			"GO:0005901", // caveola
			"GO:0005882", // intermediate filament
			"GO:0045862", // positive regulation of proteolysis [filters cln6)??
							// would remove E3
			"GO:0006614", // SRP-dependent cotranslational protein targeting to
							// membrane
			"GO:0006613", // cotranslational protein targeting to membrane
			"GO:0045047", // protein targeting to ER
			"GO:0043248", // proteasome assembly
			"GO:0017048", // Rho GTPase binding
			"GO:0035023", // regulation of Rho protein signal transduction
			"GO:0046968", // peptide antigen transport
			"GO:0005178", // integrin binding
			"GO:0098640",// integrin binding involved in cell-matrix adhesion

	};
	private final static String[] GOPartNameToExclude = { "Spliceosome", "RNA polymerase II", // "neuron",
			"DNA helicase", "Centrosome", "transcription factor activity", // "mitochondrion"
	};
	private final static String[] GOPartNameToInclude = { "myosin", "kinesin" };
	private final static GORetriever goRetriever = new GORetriever(new File("z:\\share\\Salva\\data\\go"));

	private final Set<String> filteredOut = new HashSet<String>();
	private final String name;

	@Override
	public String toString() {
		final StringBuilder sb = new StringBuilder();
		sb.append(
				"Inclusion list (Any protein annotated with any of these terms passes the filter even if having another term from the exclusion list):");
		for (final String go : GOFilter.GOToInclude) {
			final OntologyTermI goTerm = goRetriever.getGOTermByID(go);
			sb.append(goTerm.getTermAccession() + "\t" + goTerm.getPreferredName());
		}
		sb.append(
				"\nExclusion list (Any protein annotated with any of these following terms will be discarded unless having other one from the inclusion list)");
		for (final String go : GOFilter.GOToExclude) {

			final OntologyTermI goTerm = goRetriever.getGOTermByID(go);
			if (goTerm != null) {
				sb.append(goTerm.getTermAccession() + "\t" + goTerm.getPreferredName());
			}
		}
		sb.append(
				"\nExclusion list names (Any protein annotated with a GO term containing the following text will be discarded unless it is annotated with any GO term from the inclusion list)");
		for (final String name : GOFilter.GOPartNameToExclude) {

			sb.append(name);
		}
		return sb.toString();
	}

	public GOFilter(String name) {
		this.name = name;
	}

	private final Set<String> valid = new HashSet<String>();

	@Override
	public boolean isValid(String accession) {
		final String parsedAcc = FastaParser.getACC(accession).getAccession();
		if (valid.contains(parsedAcc)) {
			return true;
		} else if (filteredOut.contains(parsedAcc)) {
			return false;
		}
		filterProteinsByGO(accession);
		return isValid(accession);
	}

	public void filterProteinsByGO(String rawAccession) {
		final Set<String> accs = new HashSet<String>();
		accs.add(rawAccession);
		filterProteinsByGO(accs);
	}

	public void filterProteinsByGO(Collection<String> rawAccessions) {
		try {
			final List<String> accessions = new ArrayList<String>();
			for (final String rawAccession : rawAccessions) {
				accessions.add(FastaParser.getACC(rawAccession).getAccession());
			}

			int numDiscarded = 0;
			int numValid = 0;
			final Map<String, Set<GoEntry>> goEntries = goRetriever.retrieveGOEntries(accessions);
			final Map<String, Entry> annotatedProtein = Constants.upr.getAnnotatedUniprotEntries(accessions);
			for (final String acc : accessions) {
				if (filteredOut.contains(acc) || valid.contains(acc)) {
					continue;
				}
				boolean discard = false;
				try {
					// inclusion list
					boolean notExclude = false;
					for (final String go2 : GOToInclude) {
						if (goRetriever.containsGOTerm(acc, go2)) {
							notExclude = true;
							break;
						}

					}
					for (final String partNameToInclude : GOPartNameToInclude) {
						if (goRetriever.containsTermNamePart(acc, partNameToInclude)) {
							notExclude = true;
							break;
						}
					}
					if (notExclude) {
						continue;
					}
					// exclusion list
					for (final String go : GOToExclude) {
						if (goRetriever.containsGOTerm(acc, go)) {
							discard = true;
							break;
						}
					}
					if (discard) {
						continue;
					}

					// exclusion part name list
					for (final String partNameToExclude : GOPartNameToExclude) {
						if (goRetriever.containsTermNamePart(acc, partNameToExclude)) {
							discard = true;
							break;
						}
					}
					if (discard) {
						continue;
					}
					if (annotatedProtein.containsKey(acc)) {
						final Entry entry = annotatedProtein.get(acc);
						final List<DbReferenceType> dbReferences = entry.getDbReference();
						for (final DbReferenceType dbReference : dbReferences) {
							if (dbReference.getType().equals("GO")) {
								for (final String go : GOToExclude) {
									if (dbReference.getId().equals(go)) {
										discard = true;
										break;
									}
								}
								if (discard) {
									break;
								}
							}
						}
					}
				} finally {
					if (discard) {
						numDiscarded++;
						filteredOut.add(acc);
					} else {
						valid.add(acc);
						numValid++;
					}
					// System.out.println(acc + "\t" + !discard);
				}
			}
			if (numDiscarded != 0) {
				log.debug(filteredOut.size() + " proteins discarded (" + name + ")");
				log.debug(valid.size() + " proteins valid (" + name + ")");
				log.debug((filteredOut.size() + valid.size()) + " total proteins (" + name + ")");
			}
		} catch (final IOException e) {
			e.printStackTrace();
		}
	}

	@Override
	public boolean isValid(Protein protein) {
		return isValid(protein.getAccession());
	}

}
