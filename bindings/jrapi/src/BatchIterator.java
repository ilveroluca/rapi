
package it.crs4.rapi;

public class BatchIterator implements java.util.Iterator<Fragment> {
  private Batch batch;
  private long nextFragment;
  private long end;
  private Fragment frag;

  public BatchIterator(Batch b)
  {
    if (b == null)
      throw new NullPointerException();
    batch = b;
    nextFragment = 0;
    end = batch.getNFragments();
    frag = new Fragment(batch.getNReadsPerFrag());
  }

  public Fragment next()
  {
    if (nextFragment >= end)
      throw new java.util.NoSuchElementException();

    for (int i = 0; i < batch.getNReadsPerFrag(); ++i)
      frag.set(i, batch.getRead(nextFragment, i));

    nextFragment += 1;
    return frag;
  }

  public boolean hasNext()
  {
    return nextFragment < end;
  }

  public void remove()
  {
    throw new UnsupportedOperationException();
  }
}
