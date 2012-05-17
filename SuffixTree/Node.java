import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;


public class Node {
	
	protected CharSequence edge;
	//The string is the string that the suffix was found in (for leaves)
	//The int is the starting place of the suffix
	protected Pair<String, Integer> label;
	protected Node father;
	private boolean leftDiversity;
	private List<Node> childrenList;
	
	public Node(Pair<String, Integer> label, CharSequence edge)
	{
		this.edge = edge;
		this.label = label;
		leftDiversity = false;
		childrenList = new ArrayList<Node>();
	}
	
	public static void appendLeaf(Node father, Node child)
	{
		father.addChild(child);
	}
	
	public Pair<Integer, Node> findPath(CharSequence suffix)
	{
		boolean matchFound = false;
		int lengthFound = 0;
		CharSequence edge = getEdge();
		for (int j = 0; j < edge.length(); j++)
		{
			if (suffix.charAt(j) == edge.charAt(j))
			{
				lengthFound++;
				matchFound = true;
			}
			else
			{
				break;
			}
		}
		if (matchFound)
		{
			return new Pair<Integer, Node>(lengthFound, this);
		}
		return null;
	}
	
	public void setFather(Node father)
	{
		this.father = father;
	}
	
	public Node getFather()
	{
		return father;
	}
	
	public String getPathLabel()
	{
		String path = "";
		if (father != null)
		{
			path = father.getPathLabel();
		}
		path += getEdge();
		return path; 
	}
	
	/**
	 * Returns the left characters of all of the leaves that have this
	 * node as its ancestor.
	 * @return
	 */
	public String getLeftCharacters()
	{
		StringBuffer leftChars = new StringBuffer("");
		for (int i = 0; i < childrenList.size(); i++)
		{
			leftChars.append(childrenList.get(i).getLeftCharacters());
		}
		return new String(leftChars);
	}
	
	public List<Node> getChildren()
	{
		return childrenList;
	}
	
	public void addChild(Node node)
	{
		String leftChars = getLeftCharacters();
		String nodeLeftChars = node.getLeftCharacters();
		if (leftChars.length() != 0)
		{
			for (int i = 0; i < nodeLeftChars.length(); i++)
			{
				if (!(leftChars.contains(Character.toString(nodeLeftChars.charAt(i)))))
				{
					leftDiversity = true;
				}
			}
		}
		childrenList.add(node);
	}
	
	public void removeChild(Node node)
	{
		childrenList.remove(node);		
	}
	
	public Pair<String, Integer> getLabel()
	{
		return label;
	}
	
	public void setLabel(Pair<String, Integer> label)
	{
		this.label = label;
	}
	
	public CharSequence getEdge()
	{
		return edge;
	}
	
	public void setEdge(CharSequence edge)
	{
		this.edge = edge;
	}
	
	public boolean getLeftDiversity()
	{
		return leftDiversity;
	}
	
	public List<Node> getAncestors()
	{
		List<Node> ancestors = new ArrayList<Node>();
		ancestors.add(this);
		if (father != null)
		{
			List<Node> fatherAncestors = father.getAncestors();
			for (int i = 0; i < fatherAncestors.size(); i++)
			{
				ancestors.add(i + 1, fatherAncestors.get(i));
			}
		}
		return ancestors;
	}
}
