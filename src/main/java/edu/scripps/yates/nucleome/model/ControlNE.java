package edu.scripps.yates.nucleome.model;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import edu.scripps.yates.utilities.files.FileUtils;

public class ControlNE {
	private static Set<String> accs;

	public static Set<String> getAccs() {
		if (accs == null) {
			try {
				load();
			} catch (IOException e) {
				e.printStackTrace();
				return Collections.emptySet();
			}
		}
		return accs;
	}

	private static void load() throws IOException {
		File inputFile = new File("C:\\Users\\Salva\\Desktop\\data\\4D_Nucleome\\NE_Control_50.txt");
		accs = new HashSet<String>();
		final List<String> readColumnFromTextFile = FileUtils.readColumnFromTextFile(inputFile, "\t", 0, true);
		accs.addAll(readColumnFromTextFile);
	}

	public static boolean isControl(String longerAccession) {
		final Set<String> accs2 = getAccs();
		for (String acc : accs2) {
			if (longerAccession.contains(acc)) {
				return true;
			}
		}
		return false;
	}
}
