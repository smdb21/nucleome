package edu.scripps.yates.nucleome.turboID.annotations;

import java.io.File;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import edu.scripps.yates.annotations.uniprot.UniprotProteinRetriever;
import edu.scripps.yates.nucleome.turboID.ProteinFromTurboID;
import edu.scripps.yates.utilities.proteomicsmodel.AnnotationType;
import edu.scripps.yates.utilities.proteomicsmodel.Protein;
import edu.scripps.yates.utilities.proteomicsmodel.ProteinAnnotation;

public class AnnotationsUtil {
	private final static Logger log = Logger.getLogger(AnnotationsUtil.class);
	private final UniprotProteinRetriever upr;

	public AnnotationsUtil(UniprotProteinRetriever upr) {
		this.upr = upr;
	}

	public AnnotationsUtil(File uniprotReleasesFolder) {
		upr = new UniprotProteinRetriever(null, uniprotReleasesFolder, true);
	}

	public boolean isNucleus(ProteinFromTurboID protein) {
		return isNucleus(protein.getAcc());
	}

	public boolean isNucleus(String acc) {

		final Map<String, Protein> annotatedProteins = upr.getAnnotatedProtein(acc);

		if (annotatedProteins.containsKey(acc)) {
			final Protein protein2 = annotatedProteins.get(acc);
			if (protein2 != null) {
				final Set<ProteinAnnotation> annotations = protein2.getAnnotations();
				for (final ProteinAnnotation proteinAnnotation : annotations) {
					if (proteinAnnotation.getAnnotationType() == AnnotationType.subcellular_location) {
						if (proteinAnnotation.getName() != null
								&& proteinAnnotation.getName().toLowerCase().equals("nucleus")) {
							return true;
						}
					}
				}
				// for (final ProteinAnnotation proteinAnnotation : annotations)
				// {
				// if (proteinAnnotation.getAnnotationType() ==
				// AnnotationType.GO) {
				//
				// // try with the ontology query
				// final String parentID = "GO:0005634";// nucleus
				// final String termID = proteinAnnotation.getName();
				// final boolean valid =
				// OntologyQueryInterface.containsAsParent(termID, parentID,
				// true,
				// GO_TERM_DISTANCE);
				// if (valid) {
				// return true;
				// }
				// }
				// }
			}
		}
		return false;
	}

	public boolean isDNABinding(ProteinFromTurboID protein) {
		return isDNABinding(protein.getAcc());
	}

	public boolean isDNABinding(String acc) {

		final Map<String, Protein> annotatedProteins = upr.getAnnotatedProtein(acc);

		if (annotatedProteins.containsKey(acc)) {
			final Protein protein2 = annotatedProteins.get(acc);
			if (protein2 != null) {
				final Set<ProteinAnnotation> annotations = protein2.getAnnotations();
				for (final ProteinAnnotation proteinAnnotation : annotations) {
					if (proteinAnnotation.getAnnotationType() == AnnotationType.GO) {
						if (proteinAnnotation.getName() != null
								&& proteinAnnotation.getValue().toLowerCase().contains("dna binding")) {
							return true;
						}

					} else if (proteinAnnotation.getAnnotationType() == AnnotationType.uniprotKeyword) {
						if (proteinAnnotation.getName().toLowerCase().equals("dna-binding")) {
							return true;
						}
					}
				}

				// for (final ProteinAnnotation proteinAnnotation : annotations)
				// {
				// if (proteinAnnotation.getAnnotationType() ==
				// ) {
				//
				// // try with the ontology query
				// final String parentID = "GO:0003677";// DNA binding
				// final String termID = proteinAnnotation.getName();
				// final boolean valid =
				// OntologyQueryInterface.containsAsParent(termID, parentID,
				// true,
				// GO_TERM_DISTANCE);
				// if (valid) {
				// return true;
				// }
				// }
				// }
			}
		}
		return false;
	}

	public boolean isRNABinding(ProteinFromTurboID protein) {
		return isRNABinding(protein.getAcc());
	}

	public boolean isRNABinding(String acc) {
		final Map<String, Protein> annotatedProteins = upr.getAnnotatedProtein(acc);

		if (annotatedProteins.containsKey(acc)) {
			final Protein protein2 = annotatedProteins.get(acc);
			if (protein2 != null) {
				final Set<ProteinAnnotation> annotations = protein2.getAnnotations();
				for (final ProteinAnnotation proteinAnnotation : annotations) {
					if (proteinAnnotation.getAnnotationType() == AnnotationType.GO) {
						if (proteinAnnotation.getName() != null
								&& proteinAnnotation.getValue().toLowerCase().contains("rna binding")) {
							return true;
						}

					} else if (proteinAnnotation.getAnnotationType() == AnnotationType.uniprotKeyword) {
						if (proteinAnnotation.getName().toLowerCase().equals("rna-binding")) {
							return true;
						}
					}
				}
				// for (final ProteinAnnotation proteinAnnotation : annotations)
				// {
				// if (proteinAnnotation.getAnnotationType() ==
				// AnnotationType.GO) {
				// // try with the ontology query
				// final String parentID = "GO:0003723";// RNA binding
				// final String termID = proteinAnnotation.getName();
				// final boolean valid =
				// OntologyQueryInterface.containsAsParent(termID, parentID,
				// true,
				// GO_TERM_DISTANCE);
				// if (valid) {
				// return true;
				// }
				// }
				// }
			}
		}
		return false;
	}

	public boolean isHeterochromatin(ProteinFromTurboID protein) {
		return isHeterochromatin(protein.getAcc());
	}

	public boolean isHeterochromatin(String acc) {
		final Map<String, Protein> annotatedProteins = upr.getAnnotatedProtein(acc);

		if (annotatedProteins.containsKey(acc)) {
			final Protein protein2 = annotatedProteins.get(acc);
			if (protein2 != null) {
				final Set<ProteinAnnotation> annotations = protein2.getAnnotations();
				for (final ProteinAnnotation proteinAnnotation : annotations) {
					if (proteinAnnotation.getAnnotationType() == AnnotationType.GO) {
						if (proteinAnnotation.getName() != null
								&& proteinAnnotation.getValue().toLowerCase().contains("heterochromatin")) {
							return true;
						}

					}
				}
				// for (final ProteinAnnotation proteinAnnotation : annotations)
				// {
				// if (proteinAnnotation.getAnnotationType() ==
				// AnnotationType.GO) {
				// final String parentID = "GO:0000792";// heterochromatin
				// final String termID = proteinAnnotation.getName();
				// // try with the ontology query
				// final boolean valid =
				// OntologyQueryInterface.containsAsParent(termID, parentID,
				// true,
				// GO_TERM_DISTANCE);
				// if (valid) {
				// return true;
				// }
				// }
				// }
			}
		}
		return false;
	}

	public boolean isTranscriptionFactor(ProteinFromTurboID protein) {
		return isTranscriptionFactor(protein.getAcc());
	}

	public boolean isTranscriptionFactor(String acc) {
		final Map<String, Protein> annotatedProteins = upr.getAnnotatedProtein(acc);
		if (acc.equals("")) {
			log.info("P25976");
		}
		if (annotatedProteins.containsKey(acc)) {
			final Protein protein2 = annotatedProteins.get(acc);
			if (protein2 != null) {
				final Set<ProteinAnnotation> annotations = protein2.getAnnotations();
				for (final ProteinAnnotation proteinAnnotation : annotations) {
					if (proteinAnnotation.getAnnotationType() == AnnotationType.GO) {
						if (proteinAnnotation.getName() != null
								&& proteinAnnotation.getValue().toLowerCase().contains("transcription factor")) {
							return true;
						}

					} else if (proteinAnnotation.getAnnotationType().getKey().equals("chain")) {
						if (proteinAnnotation.getName().toLowerCase().contains("transcription factor")) {
							return true;
						}
					}

				}
			}
		}
		return false;
	}

	public String getDescription(ProteinFromTurboID protein) {
		return getDescription(protein.getAcc());
	}

	public String getDescription(String acc) {
		final Map<String, Protein> annotatedProteins = upr.getAnnotatedProtein(acc);
		if (annotatedProteins.containsKey(acc)) {
			final Protein protein2 = annotatedProteins.get(acc);
			if (protein2 != null) {
				return protein2.getDescription();
			}
		}
		return "-";
	}

	public boolean getTransmembraneRegion(ProteinFromTurboID protein) {
		return getTransmembraneRegion(protein.getAcc());
	}

	public boolean getTransmembraneRegion(String acc) {

		final Map<String, Protein> annotatedProtein = upr.getAnnotatedProtein(acc);
		if (annotatedProtein.containsKey(acc)) {
			final Protein protein2 = annotatedProtein.get(acc);
			if (protein2 != null) {
				final Set<ProteinAnnotation> annotations = protein2.getAnnotations();
				for (final ProteinAnnotation proteinAnnotation : annotations) {
					if (proteinAnnotation.getAnnotationType().getKey().equals("transmembrane region")) {
						return true;
					}
				}

			}
		}

		return false;
	}
}
