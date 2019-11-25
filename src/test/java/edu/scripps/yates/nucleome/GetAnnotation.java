package edu.scripps.yates.nucleome;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import org.junit.Test;

import edu.scripps.yates.nucleome.turboID.annotations.AnnotationsUtil;

public class GetAnnotation {

	final AnnotationsUtil annotationsUtil = new AnnotationsUtil(Constants.upr);

	@Test
	public void getTransmembraneRegion() {
		final File input = new File("Z:\\share\\Salva\\data\\4D_Nucleome\\TMT8_EMD_IP\\test.txt");
		final File output = new File("Z:\\\\share\\\\Salva\\\\data\\\\4D_Nucleome\\\\TMT8_EMD_IP\\\\test_transmem.txt");
		FileWriter fw = null;
		try {
			fw = new FileWriter(output);
			final List<String> accs = Files.readAllLines(Paths.get(input.toURI()));
			Constants.upr.getAnnotatedProteins(accs);
			for (final String acc : accs) {
				fw.write(annotationsUtil.getTransmembraneRegion(acc) + "\t" + acc + "\n");
			}
		} catch (final IOException e) {
			e.printStackTrace();
		} finally {
			try {
				fw.close();
			} catch (final IOException e) {
				e.printStackTrace();
			}
		}
	}

	@Test
	public void getAnnotations() {
		final File input = new File("Z:\\share\\Salva\\data\\4D_Nucleome\\TMT8_EMD_IP\\test.txt");
		final File output = new File(
				"Z:\\\\share\\\\Salva\\\\data\\\\4D_Nucleome\\\\TMT8_EMD_IP\\\\test_annotated.txt");

		FileWriter fw = null;
		try {
			fw = new FileWriter(output);
			final List<String> accs = Files.readAllLines(Paths.get(input.toURI()));
			Constants.upr.getAnnotatedProteins(accs);
			fw.write(
					"ACC\tTransmembrane region\tnucleus\tDNA binding\ttranscription factor\tRNA binding\theterochromatin\n");
			for (final String acc : accs) {
				fw.write(acc + "\t" + annotationsUtil.getTransmembraneRegion(acc) + "\t"
						+ annotationsUtil.isNucleus(acc) + "\t" + annotationsUtil.isDNABinding(acc) + "\t"
						+ annotationsUtil.isTranscriptionFactor(acc) + "\t" + annotationsUtil.isRNABinding(acc) + "\t"
						+ annotationsUtil.isHeterochromatin(acc) + "\n");
			}
		} catch (final IOException e) {
			e.printStackTrace();
		} finally {
			try {
				fw.close();
			} catch (final IOException e) {
				e.printStackTrace();
			}
		}
	}
}
