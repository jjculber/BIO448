
public class Leaf extends Node {
	
	private Character leftCharacter;

	public Leaf(Pair<String, Integer> label, CharSequence edge, Character leftCharacter)
	{
		super(label, edge);
		this.leftCharacter = leftCharacter;
	}
	
	@Override
	public String getLeftCharacters()
	{
		return Character.toString(leftCharacter);
	}
	
	public Character getLeftCharacter() {
		return leftCharacter;
	}

	public void setLeftCharacter(Character leftCharacter) {
		this.leftCharacter = leftCharacter;
	}
}
