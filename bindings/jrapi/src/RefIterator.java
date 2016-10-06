
package it.crs4.rapi;

import java.util.Iterator;

public class RefIterator implements Iterator<Contig> {
  private Ref ref;
  private int nextContig;
  private int nContigs;

  public RefIterator(Ref ref)
  {
    if (ref == null)
      throw new NullPointerException();
    this.ref = ref;
    nextContig = 0;
    nContigs = this.ref.getNContigs();
  }

  public Contig next()
  {
    if (nextContig >= nContigs)
      throw new java.util.NoSuchElementException();

    Contig c = ref.getContig(nextContig);
    nextContig += 1;
    return c;
  }

  public boolean hasNext()
  {
    return nextContig < nContigs;
  }

  public void remove()
  {
    throw new UnsupportedOperationException();
  }
}
