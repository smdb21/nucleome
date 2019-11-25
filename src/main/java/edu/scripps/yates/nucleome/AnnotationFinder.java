package edu.scripps.yates.nucleome;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;

import edu.scripps.yates.annotations.uniprot.UniprotProteinRetriever;
import edu.scripps.yates.utilities.proteomicsmodel.Protein;
import edu.scripps.yates.utilities.proteomicsmodel.ProteinAnnotation;

public class AnnotationFinder {
	private final File uniprotReleasesFolder = new File("Z:\\share\\Salva\\data\\uniprotKB");

	@Test
	public void annotationFinder() {
		try {
			final File file = new File(
					"D:\\Dropbox (Scripps Research)\\4DN shared folder\\TurboID\\NE4_NuCy_T11\\accs.txt");
			final FileWriter fw = new FileWriter(file.getParent() + File.separator + "annotated.txt");
			final List<String> accs = Files.readAllLines(Paths.get(file.toURI()));
			final UniprotProteinRetriever upr = new UniprotProteinRetriever(null, uniprotReleasesFolder, true);
			final Map<String, Protein> annotatedProteins = upr.getAnnotatedProteins(accs);
			for (final String acc : accs) {
				fw.write(acc + "\t");
				if (annotatedProteins.containsKey(acc)) {
					final Protein protein2 = annotatedProteins.get(acc);
					if (protein2 != null) {
						final Set<ProteinAnnotation> annotations = protein2.getAnnotations();
						boolean found = false;
						for (final ProteinAnnotation proteinAnnotation : annotations) {
							if (proteinAnnotation.getAnnotationType().getKey().equals("transmembrane region")) {
								found = true;
								break;
							}
						}
						if (found) {
							fw.write("TRUE\n");
						} else {
							fw.write("FALSE\n");
						}
					} else {
						fw.write("-\n");
					}
				} else {
					fw.write("-\n");
				}
			}
			fw.close();
		} catch (final IOException e) {
			e.printStackTrace();
		}
	}
}
