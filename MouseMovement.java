
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
//package logreader;

/**
 *
 * @author Angelica
 */
class MouseMovement extends UserEvent {
    private Double xLoc; 
    private Double yLoc;
  MouseMovement(){
      super();
      xLoc = 0.0;
      yLoc = 0.0;
  }
  
  MouseMovement(Double xLoc, Double yLoc){
	  super();
	  this.xLoc = xLoc;
	  this.yLoc = yLoc;
  }
  
  MouseMovement(int xLoc, int yLoc){
	  super();
	  this.xLoc = (double)xLoc;
	  this.yLoc = (double)yLoc;
  }

   
    /**
     * @return the xLoc
     */
    public Double getxLoc() {
        return xLoc;
    }

    /**
     * @param xLoc the xLoc to set
     */
    public void setxLoc(Double xLoc) {
        this.xLoc = xLoc;
    }

    /**
     * @return the yLoc
     */
    public Double getyLoc() {
        return yLoc;
    }

    /**
     * @param yLoc the yLoc to set
     */
    public void setyLoc(Double yLoc) {
        this.yLoc = yLoc;
    }
    
    public String toString(){
        return super.toString()+":M:"+xLoc+","+yLoc;
        
    }
}
