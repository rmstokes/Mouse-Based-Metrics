import java.util.ArrayList;


/*
 * 
 * Each Subject has it's own set of 10 Templates. Each template contains a feature set for 
 * the direction, curvature angle, and curvature distance metrics for a given instance of mouse movements.
 * 
 * 
 * 
 */

public class Subject {
	
	//instance variables
	private ArrayList <Template> template; //each subject will have 10 templates containing the features for each instance of moves
	private int subjectNumber;
	private String subjectName;
	private boolean lazyDelete;
	
	//constructor
	public Subject(int subjectNumber, String subjectName){
		this.subjectNumber = subjectNumber;
		this.subjectName = subjectName;
		template = new ArrayList<Template>();
		lazyDelete = false;
	}
	
	public void addTemplate(Template template){
		//add template of features
		//System.out.println("inside addTemplate()...");
		this.template.add(template);
	}
	
	public void addSubjectNumber(int subjectNumber){
		this.subjectNumber = subjectNumber;
	}
	
	public int getSubjectNumber(){
		
		return subjectNumber;
	}
	
	public ArrayList<Template> getTemplate(){
		
		return template;
	}
	
	public void delete(){
		lazyDelete = true;
	}
	
	public void undelete(){
		lazyDelete = false;
	}
	
	public boolean isDeleted(){
		return lazyDelete;
	}

}
