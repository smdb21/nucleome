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

import edu.scripps.yates.utilities.proteomicsmodel.Protein;
import edu.scripps.yates.utilities.proteomicsmodel.ProteinAnnotation;

public class GetAnnotation {
	private boolean getTransmembraneRegion(String acc) throws IOException {

		final Map<String, Protein> annotatedProtein = Constants.upr.getAnnotatedProtein(acc);
		if (annotatedProtein.containsKey(acc)) {
			final Protein protein2 = annotatedProtein.get(acc);
			if (protein2 != null) {
				Set<ProteinAnnotation> annotations = protein2.getAnnotations();
				for (ProteinAnnotation proteinAnnotation : annotations) {
					if (proteinAnnotation.getAnnotationType().getKey().equals("transmembrane region")) {
						return true;
					}
				}

			}
		}

		return false;
	}

	@Test
	public void getTransmembraneRegion() {
		final File input = new File(
				"Z:\\share\\Salva\\data\\4D_Nucleome\\Xi data March18\\total_protein_accs_from_Xi.txt");
		final File output = new File(
				"Z:\\share\\Salva\\data\\4D_Nucleome\\Xi data March18\\total_protein_accs_from_Xi_transmem.txt");
		FileWriter fw = null;
		try {
			fw = new FileWriter(output);
			List<String> accs = Files.readAllLines(Paths.get(input.toURI()));
			Constants.upr.getAnnotatedProteins(accs);
			for (String acc : accs) {
				fw.write(getTransmembraneRegion(acc) + "\t" + acc + "\n");
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				fw.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
}
