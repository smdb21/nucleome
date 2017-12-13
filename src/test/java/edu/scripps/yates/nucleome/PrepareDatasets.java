package edu.scripps.yates.nucleome;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

import org.junit.Test;

import edu.scripps.yates.nucleome.filters.GOFilter;
import edu.scripps.yates.utilities.files.FileUtils;
import edu.scripps.yates.utilities.remote.RemoteSSHFileReference;

public class PrepareDatasets {
	@Test
	public void testingGOFilter() {
		GOFilter filter = new GOFilter(" ");
		filter.isValid("O54962");
	}

	private final String datasetsPathsFile = "z:\\share\\Salva\\data\\4D_Nucleome\\datasets_paths.txt";
	private static final String output = "z:\\share\\Salva\\data\\4D_Nucleome\\PACom";

	@Test
	public void getDTASelects() {
		Path path = Paths.get(new File(datasetsPathsFile).toURI());
		try {
			List<String> lines = Files.readAllLines(path);
			for (String line : lines) {
				String[] split = line.split(":");
				String name = split[0];
				String remotePath = split[1].trim();
				String hostName = "jaina.scripps.edu";
				String userName = "salvador";
				String pass = "Natjeija21";
				File outputFile = new File(output + File.separator + name + ".xml");
				RemoteSSHFileReference ref = new RemoteSSHFileReference(hostName, userName, pass,
						"DTASelect-filter.txt", outputFile);
				ref.setRemotePath(remotePath);
				File outputFile2 = ref.getRemoteFile();
				if (outputFile2.exists()) {
					System.out.println(outputFile2.getAbsolutePath() + " downloaded "
							+ FileUtils.getDescriptiveSizeFromBytes(outputFile2.length()));
				} else {
					System.out.println("Error");
				}
			}

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	@Test
	public void getPACOMBatchImport() {
		File outputFile = new File(output + File.separator + "4DNucleome_batch_Import.txt");
		FileWriter fw = null;
		try {
			fw = new FileWriter(outputFile);
			int num = 1;
			for (File file : new File(output).listFiles()) {
				if (!file.isFile()) {
					continue;
				}
				if (!file.getAbsolutePath().endsWith(".xml")) {
					continue;
				}
				fw.write("START\t" + num++ + "\n");
				fw.write("PROJECT\t4DNucleome\n");
				fw.write("DTASELECT\t" + file.getAbsolutePath() + "\n");
				fw.write("END" + "\n");
			}
			fw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
