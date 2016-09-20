
package it.crs4.rapi;

import java.io.File;

public class RapiUtils
{
  public static void loadPlugin()
  {
    loadPlugin(null);
  }

  public static void loadPlugin(String path)
  {
    String soPath;
    if (path != null)
      soPath = path;
    else
      soPath = System.getProperty("jrapi.so");

    if (soPath != null)
      System.err.println("Loading shared object from cmd line argument " + soPath);
    else {
      System.err.println("jrapi.so path not provided!");
      System.err.println("Specify the property jrapi.so or pass it as a command line argument");
      System.exit(1);
    }

    // load the jrapi shared object
    String absPath = new File(soPath).getAbsolutePath();
    System.err.println("Loading library file " + absPath);
    System.load(absPath);
  }
}

// vim: set et sw=2
