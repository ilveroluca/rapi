
package it.crs4.rapi;

import java.io.File;

public class RapiUtils
{
	private static final String SHARED_OBJ_PATH = "/jrapi.so";

  public static void loadPlugin()
  {
    loadPlugin(null);
  }

	public static void loadPlugin(String path)
	{
		if (path == null)
		{
			try {
				cz.adamh.utils.NativeUtils.loadLibraryFromJar(SHARED_OBJ_PATH);
			} catch (java.io.IOException e) {
				throw new RuntimeException(e.getMessage(), e);
			}
		}
		else
		{
			loadPluginFromPath(path);
		}
	}

  public static void loadPluginFromPath(String soPath)
  {
		if (soPath == null)
			throw new IllegalArgumentException("path must not be null");

    String absPath = new File(soPath).getAbsolutePath();
		System.err.println("Loading shared object from " + absPath);
    System.load(absPath);
  }

}

// vim: set et sw=2
