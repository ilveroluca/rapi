/************************************************************************************
 * This code is published under the The MIT License.
 *
 * Copyright (c) 2016 Center for Advanced Studies,
 *                      Research and Development in Sardinia (CRS4), Pula, Italy.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ************************************************************************************/



/*
 * Run this with a command line like:
 * java -cp build/jrapi.jar -Djrapi.so=$PWD/jrapi.so rapi_example /path/to/ref.fasta < /input/data.prq > /output/path.sam
 */

import it.crs4.rapi.AlignerState;
import it.crs4.rapi.Batch;
import it.crs4.rapi.Opts;
import it.crs4.rapi.Rapi;
import it.crs4.rapi.RapiConstants;
import it.crs4.rapi.RapiException;
import it.crs4.rapi.RapiUtils;
import it.crs4.rapi.Ref;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.IOException;

public class rapi_example
{
  private Opts opts;
  private String refPath;
  private AlignerState aligner;
  private int linesRead = 0;
  private SimpleLogger log = new SimpleLogger();

  private static final int BATCH_SIZE = 5000;

  public rapi_example() throws RapiException
  {
    RapiUtils.loadPlugin();
    opts = new Opts();
    opts.setShareRefMem(true);
    opts.setNThreads(8);
    Rapi.init(opts);
    aligner = new AlignerState(opts);
  }

  public class SimpleLogger
  {
    private int loglevel = 5;

    public SimpleLogger(int level) {
      loglevel = level;
    }

    public SimpleLogger() {
      this(5);
    }

    protected void log(int level, String fmt, Object[] args)
    {
      try {
        if (level >= loglevel) {
          System.err.format(fmt, args);
          System.err.println();
        }
      }
      catch (Exception e) {
        System.err.println("Exception raised while trying to log: " + e);
        System.err.println("fmt: " + fmt);
      }
    }

    public void debug(String fmt, Object... args) {
      log(5, "DEBUG: " + fmt, args);
    }
  }

  protected void close() throws RapiException
  {
    Rapi.shutdown();
  }

  protected boolean loadBatch(BufferedReader in, Batch dest) throws RapiException, IOException
  {
    if (BATCH_SIZE <= 0)
      throw new RuntimeException("BUG!  BATCH_SIZE must be > 0");
    dest.reserve(BATCH_SIZE);
    dest.clear();

    int nLines = 0;
    String thisLine = null;

    while (nLines < BATCH_SIZE)
    {
      thisLine = in.readLine();
      if (thisLine == null) {
        log.debug("input EOF");
        break;
      }

      nLines += 1;
      String[] parts = thisLine.split("\t");
      if (parts.length != 5)
        throw new IllegalArgumentException("Invalid prq format in line " + (nLines + linesRead) + ":\n" + thisLine);
      dest.append(parts[0], parts[1], parts[2], RapiConstants.QENC_SANGER);
      dest.append(parts[0], parts[3], parts[4], RapiConstants.QENC_SANGER);
    }

    linesRead += nLines;
    log.debug("Added %d lines to batch", nLines);
    return nLines > 0;
  }

  protected void processAlignments(Batch reads) throws RapiException
  {
    System.out.append(Rapi.formatSamBatch(reads));
  }


  protected void parseArgs(String[] args)
  {
    if (args.length != 1)
      throw new RuntimeException("Usage: cat prq | <Program> REFERENCE");
    refPath = args[0];
  }

  public void run(String[] args) throws RapiException, IOException
  {
    parseArgs(args);
    Batch reads = new Batch(2);
    log.debug("Created batch");
    Ref ref = new Ref(refPath);
    log.debug("Loaded reference from " + refPath);

    BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));

    log.debug("Starting to process");
    long startTime = System.nanoTime();

    System.out.append(Rapi.formatSamHdr(ref));

    while (loadBatch(reader, reads)) {
      log.debug("Loaded batch.  Going to align");
      aligner.alignReads(ref, reads);
      log.debug("Alignment finished.  Need to print output");
      processAlignments(reads);
      log.debug("Processed %d lines.  Time so far: %.3f s", linesRead, (System.nanoTime() - startTime) / 1.0e9);
    }

    long endTime = System.nanoTime();
    log.debug("Processed %d lines (therefore %d reads)", linesRead, linesRead * 2);
    log.debug("Total time: %.3f s", (endTime - startTime) / 1.0e9);
  }

  public static void main(String[] args) throws RapiException, IOException
  {
    rapi_example instance = new rapi_example();
    try {
      instance.run(args);
    } finally {
      instance.close();
    }
  }
}
