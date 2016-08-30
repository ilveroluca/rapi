
import it.crs4.rapi.lowrapi.*;

import java.io.File;

import org.junit.*;
import static org.junit.Assert.*;

public class TestLowRapi
{
  private static final String RELATIVE_MINI_REF = "../../tests/mini_ref/mini_ref.fasta";
  private static File jrapiBaseDir;

  private opts rapiOpts;
  private int error;
  private File miniRefPath;

  @BeforeClass
  public static void initSharedObj()
  {
    String basedirProperty = System.getProperty("jrapi.basedir");
    if (basedirProperty == null)
    {
      System.err.println("jrapi.basedir property is not defined.  Using CWD as base directory");
      jrapiBaseDir = new File(".");
    }
    else
    {
      System.err.println("Using " + basedirProperty + " as base directory");
      jrapiBaseDir = new File(basedirProperty);
    }

    loadPlugin(null);
  }

  @Before
  public void init()
  {
    miniRefPath = new File(jrapiBaseDir, RELATIVE_MINI_REF);
    rapiOpts = new opts();
    error = Rapi.init(rapiOpts);
    assertEquals(Rapi.NO_ERROR, error);
  }

  @After
  public void tearDown()
  {
    error = Rapi.shutdown();
    assertEquals(Rapi.NO_ERROR, error);

  }

  @Test
  public void testVersions()
  {
    assertTrue(Rapi.aligner_name().startsWith("bwa-mem"));
    assertTrue(Rapi.aligner_version().startsWith("0."));
    assertTrue(Rapi.plugin_version().startsWith("0."));
  }

  @Test
  public void testLoadUnloadRef()
  {
    ref refObj = new ref();
    error = Rapi.ref_load(miniRefPath.getAbsolutePath(), refObj);
    assertEquals(Rapi.NO_ERROR, error);

    assertEquals(1, refObj.getN_contigs());

    error = Rapi.ref_free(refObj);
    assertEquals(Rapi.NO_ERROR, error);
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
