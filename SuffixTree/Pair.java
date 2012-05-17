public class Pair<F,S> {

  public final F first;
  public final S second;

  public Pair(F first, S second) {
    this.first = first;
    this.second = second;
  }

  @Override
  public int hashCode() { return first.hashCode() ^ second.hashCode(); }

  @Override
  public boolean equals(Object o) {
    if (o == null) return false;
    if (!(o instanceof Pair)) return false;
    Pair pairo = (Pair) o;
    return this.first.equals(pairo.first) &&
           this.second.equals(pairo.second);
  }

}