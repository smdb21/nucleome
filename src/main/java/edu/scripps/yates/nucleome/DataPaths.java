package edu.scripps.yates.nucleome;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import edu.scripps.yates.utilities.util.Pair;

public class DataPaths {
	private final File dataPathsFile;
	private boolean loaded = false;
	private final Map<String, String> paths = new HashMap<String, String>();

	public DataPaths(String dataPathsFileName) throws FileNotFoundException {
		dataPathsFile = new File(dataPathsFileName);
		if (!dataPathsFile.exists()) {
			throw new FileNotFoundException("File " + dataPathsFileName + " is not found");
		}
	}

	public String[] getFiles(String key) throws IOException {
		if (!loaded) {
			load();
		}
		int numreplicate = -1;
		if (key.length() == 3) {
			numreplicate = Integer.valueOf(key.substring(2, 3));
			key = key.substring(0, 2);

		}
		final List<String> ret = new ArrayList<String>();
		for (final String key2 : paths.keySet()) {
			if (numreplicate == -1 && key2.toLowerCase().startsWith(key.toLowerCase())) {
				ret.add(key2);
			} else if (numreplicate != -1 && key2.toLowerCase().startsWith(key.toLowerCase())
					&& key2.endsWith(String.valueOf(numreplicate))) {
				ret.add(key2);
			}

		}
		Collections.sort(ret, new Comparator<String>() {

			@Override
			public int compare(String o1, String o2) {
				// N, Ne, C
				int num1 = 0;
				int num2 = 0;
				if (o1.contains("Ne")) {
					num1 = 1;
				} else if (o1.contains("C")) {
					num1 = 2;
				}
				if (o2.contains("Ne")) {
					num2 = 1;
				} else if (o2.contains("C")) {
					num2 = 2;
				}
				return Integer.compare(num1, num2);
			}
		});
		final String[] array = ret.toArray(new String[0]);
		final String[] arrayPaths = new String[array.length];
		for (int i = 0; i < array.length; i++) {
			arrayPaths[i] = paths.get(array[i]);
		}
		return arrayPaths;
	}

	public Pair<String, String> getXiFiles(String key) throws IOException {
		if (!loaded) {
			load();
		}

		for (final String key2 : paths.keySet()) {
			if (key2.contains(key)) {
				final String string = paths.get(key2);
				return new Pair<String, String>(key, string);
				// return string.replace("/data/2/rpark/ip2_data//",
				// "/ip2_garibaldi/");
			}
		}
		return null;
	}

	private void load() throws IOException {
		final Stream<String> linesStream = Files.lines(dataPathsFile.toPath(), Charset.defaultCharset());
		final List<String> lines = linesStream.collect(Collectors.toList());
		for (String line : lines) {
			line = line.trim();
			if ("".equals(line)) {
				continue;
			}
			if (line.contains(":")) {
				final String[] split = line.split("\t");
				final String key = split[2].trim();
				if (split.length > 6) {
					final String path = split[6].trim();
					if (!"".equals(key)) {
						paths.put(key, path);
					}
				}
			}
		}
		linesStream.close();
		loaded = true;
	}
}
