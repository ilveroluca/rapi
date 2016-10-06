
package it.crs4.rapi;

import java.util.Iterator;

public class ReadAlnIterator implements Iterator<Alignment> {
  private Read read;
  private int nextAln;
  private int nAlns;

  public ReadAlnIterator(Read r)
  {
    if (r == null)
      throw new NullPointerException();
    read = r;
    nextAln = 0;
    nAlns = read.getNAlignments();
  }

  public Alignment next()
  {
    if (nextAln >= nAlns)
      throw new java.util.NoSuchElementException();

    Alignment a = read.getAln(nextAln);
    nextAln += 1;
    return a;
  }

  public boolean hasNext()
  {
    return nextAln < nAlns;
  }

  public void remove()
  {
    throw new UnsupportedOperationException();
  }
}
