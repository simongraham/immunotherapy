'''Simon Graham Immunology Model

Model is initiated with an initial population of T-cells and cancer cells. This diversiy can be regulated by
changing sigma for both cancer and T-cell population when sampling from the log-normal distribution. 

Here, there is no adoptive process and all sites are filled up correponding to our initial sampling from the log-normal distribution.

T-cells move by swapping places with other T-cells and T-cells kill cancer cells corresponding to a given exponential function.
'''

from __future__ import division
import Tkinter
import numpy 
import random
import math


size =   60       
squareWidth = 9           
BLACK = (0,0,0)
BLUE = (0,0,225)
RED = (255,0,0)
GREEN = (0,255,0)
BROWN = (165,42, 42)
OLIVE = (107, 142, 35)
t = 0

numClonotypes = 40 
numCancerCells = 40
numHealthyCells = 40
numTCells = 40
numDiffCancerAntigens = numCancerCells
#threshold = 0.7

''''''''''''''''''''''''''''''''''''''''''''''''
#Variables to change

sigma_tCells = 0.01
sigma_cancer = 0.01
Lambda = 50    # Cross reactivity constant
''''''''''''''''''''''''''''''''''''''''''''''''

canvasWidth = size * squareWidth    # full width of canvas in pixels

sampleSize = 10    # Won't need to take a sample here, as not using adpotive immunotherapy
freeSpaces = (size**2 - numTCells - numHealthyCells)
growthFactor = freeSpaces / sampleSize

s = [ [ [0, 0, 0, 0] for row in range(size)] for column in range(size) ]


directions = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1), (1, -1), (-1, 1)]
directions2 = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1), (1, -1), (-1, 1), (0,0)]
        
running = False

theWindow = Tkinter.Tk()            # create the GUI window
theWindow.title("Immunology Model")
theWindow.geometry('+50+50')        # get the window away from the corner


# Here's the Canvas where we draw the lattice using a Tkinter PhotoImage:
theCanvas = Tkinter.Canvas(theWindow, width=canvasWidth, height=canvasWidth)
theCanvas.pack()                    # put it at the top of the window
theImage = Tkinter.PhotoImage(width=canvasWidth, height=canvasWidth)
theCanvas.create_image((3, 3), image=theImage, anchor="nw", state="normal")
# The coordinates (3, 3) are a kludge to eliminate a mysterious offset that occurs otherwise.



# Function called when Start/Stop button is pressed:
def startStop():
    global running
    running = not running
    if running:
        goButton.config(text="Pause")
    else:
        goButton.config(text="Resume")




# Create the GUI controls:
controlFrame = Tkinter.Frame(theWindow)        # a frame to hold the GUI controls
controlFrame.pack()                            # put it below the canvas
tLabel = Tkinter.Label(controlFrame, text="Cancer Growth Rate: ")
tLabel.pack(side="left")
tSlider = Tkinter.Scale(controlFrame, from_=0, to=1, resolution=0.01, length=180, orient="horizontal")
tSlider.pack(side="top")
tSlider.set(0.2)                              
spacer = Tkinter.Frame(controlFrame, width=100)
spacer.pack(side="left")
goButton = Tkinter.Button(controlFrame, text="Start", width=12, command=startStop)
goButton.pack(side="left")
		


def initialise():

    # Sample number of t cells from each clonotype using the lognormal distribution

    mean_tCell = 1
    mu_tCell = mean_tCell
    mean_cancer = 1
    mu_cancer = mean_cancer

    v_tCell = numpy.random.lognormal(1, sigma_tCells, numClonotypes)   


    for i in range(numClonotypes):
        v_tCell[i] = math.log1p(v_tCell[i])
        v_tCell[i] = round(v_tCell[i])

    growth_tcell = numTCells/numClonotypes
    for i in range(numClonotypes):
        v_tCell[i] = round((v_tCell[i] * growth_tcell))


    diff = int(sum(v_tCell) - numTCells)

    if diff>0:
        for i in range(diff):
            n = 0
            while n == 0:
                x = random.randint(0,numClonotypes-1)
                if v_tCell[x] != 0:
                    v_tCell[x] -= 1
                    n = 1



    # Sample the number of cancer cells with a specific antigen by sampling from the lognormal distributon.

    v_cancer = numpy.random.lognormal(mu_cancer, sigma_cancer, numDiffCancerAntigens)   


    for i in range(numDiffCancerAntigens):
        v_cancer[i] = math.log1p(v_cancer[i])
        v_cancer[i] = round(v_cancer[i])

    growth_cancer = numCancerCells/numDiffCancerAntigens
    for i in range(numDiffCancerAntigens):
        v_cancer[i] = round((v_cancer[i] * growth_cancer))


    diff = int(sum(v_cancer) - numCancerCells)

    if diff>0:
        for i in range(diff):
            n = 0
            while n == 0:
                x = random.randint(0,numDiffCancerAntigens-1)
                if v_cancer[x] != 0:
                    v_cancer[x] -= 1
                    n = 1

    # Now we populate the initial grid with cancer cells and T cells with given diversity

    for tCell in range(numClonotypes):       
        for item in range(int(v_tCell[tCell])):
            n = 0
            while n == 0:
                i = random.randint(0,size-1)
                j = random.randint(0,size-1)
                if s[i][j][0] == 0:
                    s[i][j][0] = -1
                    s[i][j][2] = tCell + 1
                    n = 1

    
    for cancerCell in range(numDiffCancerAntigens):      
        for item in range(int(v_cancer[cancerCell])):
            n = 0
            while n == 0:
                i = random.randint(0,size-1)
                j = random.randint(0,size-1)
                if s[i][j][0] == 0 and s[i][j][1] != 1:
                    s[i][j][1] = 1
                    s[i][j][3] = cancerCell + 1
                    n = 1

    for healthyCell in range(numHealthyCells):
        n = 0
        while n == 0:
            i = random.randint(0,size-1)
            j = random.randint(0,size-1)
            if s[i][j][0] == 0 and s[i][j][1] == 0:
                s[i][j][1] = 2
                n=1


    #Repopulating the grid corresponding to a taken sample

	# sample = []
	# for sampleCell in range(sampleSize):
	# 	n = 0
	# 	while n == 0:
	# 		i = random.randint(0, size-1)
	# 		j = random.randint(0, size-1)
	# 		if s[i][j][0] == -1:
	# 			sample.append(s[i][j][2])
	# 			n=1


 #    for i in range(size):
 #        for j in range(size):
 #            if s[i][j][0] == 0 and s[i][j][1] == 0:
 #                element = (random.randint(1,sampleSize))-1
 #                s[i][j][0] = -1
 #                s[i][j][2] = sample[element]


    for i in range(size):
        for j in range(size):
            colorSquare(i,j)

	


# Function to color the square representing site (i,j):
def colorSquare(i, j):
    
    if s[i][j][0] == 0 and s[i][j][1] == 1:
        theColor = "RED"
    elif s[i][j][0] == -1 and s[i][j][1] == 0:
        theColor = "BLUE"
    elif s[i][j][0] == 0 and s[i][j][1] == 0:
        theColor = "BLACK"
    elif s[i][j][0] == 0 and s[i][j][1] == 2:
        theColor = "GREEN"
    elif s[i][j][0] == -1 and s[i][j][1] == 1:
    	theColor = "BROWN"
    elif s[i][j][0] == -1 and s[i][j][1] == 2:
    	theColor = "OLIVE"
    theImage.put(theColor, to=(i*squareWidth,j*squareWidth,(i+1)*squareWidth,(j+1)*squareWidth))



def isIn(i, j):
    return 0 <= i and i < size and \
            0 <= j and j < size



def diffBetweenReceptors(i, j, new_i, new_j):
    if numClonotypes and numDiffCancerAntigens == 1:
        diff = abs(s[i][j][2] - s[new_i][new_j][3])
        sim = 1-diff
        cross_reactivity = sim #*threshold
    else:
        x = (s[i][j][2]) - 1
        y = (s[new_i][new_j][3]) - 1
        z = (max(numClonotypes, numDiffCancerAntigens)) - 1
        sim = (abs(x-y))/z
        cross_reactivity = (math.exp(-Lambda*sim)) #*threshold   
        
    return cross_reactivity



def tryToMoveTCell(i, j):
        assert isIn(i, j)
        assert s[i][j][0] == -1
        
        # if sumFreeNeighbours(i,j) > 0:
        #     while True:
        k1 = random.randint(0,8)

        new_i = i + directions2[k1][0]
        new_j = j + directions2[k1][1]

        if isIn(new_i, new_j):
            if s[new_i][new_j][0] == -1:
                s[new_i][new_j][2] = s[i][j][2]
                s[i][j][2] = s[new_i][new_j][2]
                colorSquare(new_i, new_j)
                colorSquare(i, j)

            if s[new_i][new_j][0] == 0:
            	s[new_i][new_j][2] = s[i][j][2]
            	s[new_i][new_j][0] = -1
            	s[i][j][0] = 0
            	s[i][j][2] = 0
            	colorSquare(new_i, new_j)
                colorSquare(i, j)

                

            
def tryToGrowCancerCell(i, j):
        assert isIn(i, j)
        assert s[i][j][1] == +1

    	k1 = random.randint(0,7)

        new_i = i + directions[k1][0]
        new_j = j + directions[k1][1]

        if isIn(new_i, new_j):
            if s[new_i][new_j][1] != 1 and s[new_i][new_j][1] != 2:
                s[new_i][new_j][3] = s[i][j][3]     # No mutations
                s[new_i][new_j][1] = 1
                colorSquare(new_i, new_j)


def tryToGrowHealthyCell(i, j):
    assert isIn(i, j)
    assert s[i][j][1] == +2

    k1 = random.randint(0,7)

    new_i = i + directions[k1][0]
    new_j = j + directions[k1][1]

    if isIn(new_i, new_j):
        if s[new_i][new_j][1] != 1 and s[new_i][new_j][1] != 2:   
            s[new_i][new_j][1] = 2     

                

def getRandomCellByType1(cell_type):
    quantity = 0
    for i in range(0, size):
        for j in range(0, size):
            if s[i][j][1] == cell_type:
                quantity += 1

    assert quantity > 0

    selected_id = random.randint(1, quantity)

    quantity = 0
    for i in range(0, size):
        for j in range(0, size):
            if s[i][j][1] == cell_type:
                quantity += 1
                if quantity == selected_id:
                    return (i, j)
    assert False
                                    

def moveTCells():
    for i in range(size):
        for j in range(size):
            if s[i][j][0] == -1:
                tryToMoveTCell(i, j)
                

def growRandomCancerCell():
    cell = getRandomCellByType1(+1)
    tryToGrowCancerCell(cell[0], cell[1])

def growRandomHealthyCell():
    cell = getRandomCellByType1(+2)
    if cell != 0:
        tryToGrowHealthyCell(cell[0], cell[1])


def killCancerCell():
    for i in range(size):
        for j in range(size):
            if s[i][j][0] == -1:
                for d in directions2:
                    new_i = i + d[0]
                    new_j = j + d[1]

                    if isIn(new_i, new_j) and s[new_i][new_j][1] == 1:
                        # Takes into account threhsold level and cross-reactivity
                        x = diffBetweenReceptors(i, j, new_i, new_j)
                        if random.random() < x:

                            s[new_i][new_j][1] = 0
                            s[new_i][new_j][3] = 0
                            colorSquare(new_i, new_j)




def killHealthyCell():
    for i in range(size):
        for j in range(size):
            if s[i][j][0] == -1:
                for d in directions2:
                    new_i = i + d[0]
                    new_j = j + d[1]

                    if isIn(new_i, new_j) and s[new_i][new_j][1] == 2:
                        
                        # Takes into account threhsold level and cross-reactivity
                        # Cancer kill is 10 times as likely
                        if random.random() < math.exp(-Lambda*(0.01*numClonotypes)):
                            s[new_i][new_j][1] = 0
                            s[new_i][new_j][3] = 0
                            colorSquare(new_i, new_j)



# Diversity of initial T-Cell population. 0 is low diversity, 1 is high diversity

divArray = numpy.zeros(numClonotypes)
def diversity_tCell():
    for i in range(size):
        for j in range(size):
            if s[i][j][0] == -1:
                divArray[(s[i][j][2])-1] += 1

    y = len(divArray)
    counter = 0

    for i in range(y):
        z = divArray[i] * (divArray[i] - 1)
        counter += z 

    simpsons = counter / ((freeSpaces + numTCells) * ((freeSpaces + numTCells) - 1))

    return 1 - simpsons


# Diversity of initial cancer cell population. 0 is low diversity, 1 is high diversity

divArray2 = numpy.zeros(numDiffCancerAntigens)
def diversity_cancer():
    for i in range(size):
        for j in range(size):
            if s[i][j][1] == 1:
                divArray2[(s[i][j][3])-1] += 1

    y = len(divArray2)
    counter = 0

    for i in range(y):
        z = divArray2[i] * (divArray2[i] - 1)
        counter += z 

    simpsons = counter / ((freeSpaces + numCancerCells) * ((freeSpaces + numCancerCells) - 1))

    return 1 - simpsons
            
def healthyGrowthRate(numberHealthyCells):
    growthrate = 0
    if numberHealthyCells == 40:
        growthrate = 0
    elif numberHealthyCells > 30:
        growthrate = cancerGrowthProb/4
    elif numberHealthyCells > 20:
        growthrate = cancerGrowthProb/2
    elif numberHealthyCells > 10:
        growthrate = (3*cancerGrowthProb)/4
    elif numberHealthyCells > 0: 
        growthrate = cancerGrowthProb

    return growthrate

def countCancerCells():
    counter = 0
    for i in range(size):
        for j in range(size):
            if s[i][j][0] == 1:
                counter += 1
    return counter     



def countHealthyCells():
    counter = 0
    for i in range(size):
        for j in range(size):
            if s[i][j][1] == 2:
                counter += 1
    return counter

        
def simulate():
    global t
    global running
    global cancerGrowthProb
    if running:
        t += 1
        cancerGrowthProb = tSlider.get() # Set cancer growth rate
        for step in range(1):
            x = countHealthyCells()
            healthyGrowthProb = healthyGrowthRate(x)
            random.seed(t)
            prob_cancer = random.random()
            prob_healthy = random.random()


            if prob_cancer < cancerGrowthProb:
                growRandomCancerCell()
            else:
                moveTCells()

            if prob_healthy < healthyGrowthProb:
                growRandomHealthyCell()

        killCancerCell()
        killHealthyCell()
                                                 
    theWindow.after(1,simulate)       



initialise()
simulate()
theWindow.mainloop() 
