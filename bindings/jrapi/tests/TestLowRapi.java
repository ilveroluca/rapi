
import it.crs4.rapi.lowrapi.*;

import java.io.File;

import org.junit.*;
import static org.junit.Assert.*;

public class TestLowRapi
{
  private opts rapiOpts;
  private int error;

  @BeforeClass
  public static void initSharedObj()
  {
    loadPlugin(null);
  }

  @Before
  public void init()
  {
    rapiOpts = new opts();
    error = Rapi.init(rapiOpts);
    assertEquals(0, error);
  }

  @Test
  public void testVersions()
  {
    assertTrue(Rapi.aligner_name().startsWith("bwa-mem"));
    assertTrue(Rapi.aligner_version().startsWith("0."));
    assertTrue(Rapi.plugin_version().startsWith("0."));
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


  public static void main(String args[])
  {
    if (args.length == 1)
      loadPlugin(args[0]);
    else
      loadPlugin(null);

    org.junit.runner.JUnitCore.main(TestLowRapi.class.getName());
  }
}
