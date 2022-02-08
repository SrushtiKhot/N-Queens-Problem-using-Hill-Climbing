# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 08:44:33 2021

@author: srevadig
"""

import random
import numpy as np
import copy

class hillclimb:
    
    def __init__(self, state = None, m_sideways = 0, nqueens = 0):
        self.first_state = state
        
        
        if(state == None and nqueens):
            self.nqueens = nqueens
        else:
            l_s=len(state)
            self.nqueens =l_s


        self.total_number_steps = 0
        self.m_sideways = m_sideways
        self.r_sides = m_sideways
    
    
    # retrive_right_diagonal -  Outputs the diagonal cells that are to the right of the current cell.    
    def retrive_right_diagonal(self, r, c):
        i = c+1
        right_di_cells = []
        while i < self.nqueens:
            c1=r-(i-c)
            if c1 >= 0: 
                right_di_cells.append((c1,i))
            c2= r+(i-c)  
            if c2 <= self.nqueens-1: 
                right_di_cells.append((c2, i))
            i=i+1
        return right_di_cells
    
    
    # retrieve_right_horizontal - Outputs the horizontal cells that are to the right of the current cell 

    def retrieve_right_horizontal(self, r, c):
        i = c+1    
        right_h_cells = []
        while i < self.nqueens:
            right_h_cells.append((r, i))
            i+=1
        return right_h_cells

    
    #combine - Combines right horizontal cells and right diagonal cells to get all cells that are to the right of current cell.
    
    def combine(self, r, c):
        rrh=self.retrieve_right_horizontal(r,c)
        rrd=self.retrive_right_diagonal(r,c)
        return_ele=rrh+rrd
        return return_ele
    
    # queen_cell - retrieve the position of the cell of the queen of the specified tree.

    def queen_cell(self, state):
        q_cellls = []
        for columns, r in enumerate(state):
            q_cellls.append((r,columns))
        return q_cellls


    # heuristic_n_queens - Evaluate  the heuristic value for specified state.
    
    def heuristic_n_queens(self, cell_of_state):
        hv = 0
        for row,column in cell_of_state:
            x_coord = set(cell_of_state)
            y_coord = set(self.combine(row,column))
            coord_intersection = x_coord.intersection(y_coord)
            l_i=len(coord_intersection)
            hv =hv+ l_i
        return hv


    # print_state - Print the problem as matrix.
    
    def print_state(self, cell_of_state):
        
        print(cell_of_state)
        for i in range(self.nqueens):
            mat = '|'
            for j in range(self.nqueens):
                if((i,j) in cell_of_state):
                    mat =mat+ 'Q|'
                else:
                    mat =mat+ '-|'
            print(mat)

    '''
    cal_heuristic - 
        Calculate heuristic values for all the right to take the next step.
        Returns ['matrix of heuristics value','least value of Heuristics','two arrays containing row and col of right cells containing lower value of Heuristics' ]
    '''
    def cal_heuristic(self, cell_of_state):
        mzero= np.zeros((self.nqueens,self.nqueens), int)
        matix_of_heuristic =mzero + -1
        least_of_heiristic = sum(range(self.nqueens)) + 1

        for (x_coord,y_coord) in cell_of_state:
            for i in range(self.nqueens):
                if(x_coord == i):
                    pass
                else:
                    next_state = copy.deepcopy(cell_of_state)
                    tup1=(i,y_coord)
                    next_state[y_coord] =tup1
                    mh=self.heuristic_n_queens(next_state)
                    matix_of_heuristic[i,y_coord] = mh
                    
                    min_val= min(least_of_heiristic, matix_of_heuristic[i,y_coord])

                    least_of_heiristic =min_val
                    ind=np.where(matix_of_heuristic == least_of_heiristic)

        return matix_of_heuristic, least_of_heiristic,ind

    '''
    steepest_ascent_algo -Calculates the least value of heuristic and executes.we get 1 if there is Flat,shoulder 
    or flat local maxima. 2-if local maxima occurred. 3- if success occurs.
    '''
    
    def steepest_ascent_algo(self, state = None, hv = None, steps = 0):
        
        cell_of_state = None
        
        if(steps == 0):
            state = self.first_state
            cell_of_state = self.queen_cell(state)
            hv = self.heuristic_n_queens(cell_of_state)
        else:
            cell_of_state = self.queen_cell(state)
        
        steps=steps+1
        self.total_number_steps+=1

        if(hv == 0):
            print("***Success achieved***")
            self.print_state(cell_of_state)
            return 3, steps
        
        if(steps == 1):
            print("***Initial representation***")
            self.print_state(cell_of_state)
        else:
            print('Step number-->', steps)
            self.print_state(cell_of_state)
            
        matix_of_heuristic = self.cal_heuristic(cell_of_state)
        least_of_heiristic = matix_of_heuristic[1]
        
        l_mh=len(matix_of_heuristic[2][0])
        
        random_int = random.randint(0, l_mh-1)
        row = matix_of_heuristic[2][0][random_int]
        column = matix_of_heuristic[2][1][random_int]

        next_state = copy.deepcopy(state)
        next_state[column] = row

        if(least_of_heiristic < hv):
            return self.steepest_ascent_algo(next_state, least_of_heiristic, steps)
        elif (least_of_heiristic > hv): 
            print("Search procedure Failed")
            return 2, steps 
        elif (least_of_heiristic == hv): 
            if(self.r_sides): 
                self.r_sides=self.r_sides-1
                return self.steepest_ascent_algo(next_state, least_of_heiristic, steps)
            else:
                print("Search procedure Failed")
                return 1, steps

    # random_state - new random state is generated whenever this method  is called.
    
    def random_state(self):
        rset = []
        for i in range(self.nqueens):
            r_int=random.randint(0,self.nqueens-1)
            rset.append(r_int)
        return rset
        
    
    # random_restart_hc - Using Steepest Ascent as the base ,Random restart hill climbing search is implemented.
        
    def random_restart_hc(self):
        r_val = 0
        while True:
            r_val=r_val+1
            self.first_state = self.random_state()
            output = self.steepest_ascent_algo()
            if(output[0] == 3):
                t_steps=self.total_number_steps
                return r_val, output[1], t_steps
                break
    
class Hill_Climbing_check:
    
    def __init__(self, n_value, max_iter, m_sideways = 0):
        self.n_value = n_value
        self.max_iter = max_iter
        self.m_sideways = m_sideways
        self.steepest_ascent_config = [[0,[]],[0,[]],[0,[]],[0,[]]]
        self.steepest_ascent_side_config = [[0,[]],[0,[]],[0,[]],[0,[]]]
        self.random_restart_config = [0, [], [], []]
        self.random_restart_side_config = [0, [], [], []]
    
    
    # explore - Executes steepest Ascent algo and random restart hill climbing for max number of iterations times
    
    def explore(self):
        
        if(self.n_value in range(4)):
            print('VALUE MUST BE GREATER THAN 3.')
            return
        
        if(self.max_iter < 1):
            print('VALUE MUST BE GREATER THAN 1.')
            return

        for n in range(self.max_iter):
            self.steepest_ascent_config[0][0]=self.steepest_ascent_config[0][0]+1
            self.steepest_ascent_side_config[0][0]=self.steepest_ascent_side_config[0][0]+1
            self.random_restart_config[0]=self.random_restart_config[0]+1
            self.random_restart_side_config[0]= self.random_restart_side_config[0]+1
            s = []
            
            for i in range(self.n_value):
                rand_int1=random.randint(0,self.n_value-1)
                s.append(rand_int1)

            print("HILL CLIMBING SEARCH")
            hillClimbing = hillclimb(s)
            outcome = hillClimbing.steepest_ascent_algo()
            self.steepest_ascent_config[outcome[0]][0]=self.steepest_ascent_config[outcome[0]][0]+1 
            self.steepest_ascent_config[outcome[0]][1].append(outcome[1]) 
            
            print("HILL CLIMBING WITH SIDEWAYS")
            hillClimbing = hillclimb(s, self.m_sideways)
            outcome = hillClimbing.steepest_ascent_algo()
            self.steepest_ascent_side_config[outcome[0]][0]=   self.steepest_ascent_side_config[outcome[0]][0]+1 
            self.steepest_ascent_side_config[outcome[0]][1].append(outcome[1]) 
            
            print("RANDOM RESTART HILLCLIMBING")
            hillClimbing = hillclimb(None, 0, self.n_value)
            outcome = hillClimbing.random_restart_hc()
            self.random_restart_config[1].append(outcome[0]) 
            self.random_restart_config[2].append(outcome[1]) 
            self.random_restart_config[3].append(outcome[2]) 

            print("RANDOM RESTART WITH SIDEWAYS")
            hillClimbing = hillclimb(None, self.m_sideways, self.n_value)
            outcome = hillClimbing.random_restart_hc()
            self.random_restart_side_config[1].append(outcome[0]) 
            self.random_restart_side_config[2].append(outcome[1]) 
            self.random_restart_side_config[3].append(outcome[2])
        
        self.final_outcomes()
        
   
    # final_outcomes -End analysis of all 4 algorithms are printed.
    
    def final_outcomes(self):
        self.display_steepest_ascent_config(self.steepest_ascent_config, "HILL CLIMBING SEARCH")
        self.display_steepest_ascent_config(self.steepest_ascent_side_config, "HILL CLIMBING WITH SIDEWAYS")
        self.display_random_restart_config(self.random_restart_config, "RANDOM RESTART HILLCLIMBING")
        self.display_random_restart_config(self.random_restart_side_config, "RANDOM RESTART WITH SIDEWAYS")

    
    #display_random_restart_config(self) - Prints the Random Restart hill climbing search inspection.

    def display_random_restart_config(self, result, first):
        
        sum_runs = result[0]
        s1=sum(result[1])
        average_of_restart = s1 / sum_runs
        s2= sum(result[2])
        avg_last_steps = s2/ sum_runs
        s3=sum(result[3])
        avg_total = s3 / sum_runs
        
        print("\n\n"+first)
        border = ''
        lenf=len(first)
        for i in range(lenf):
            border+="-"
        print(border)
        print()
        print("n is ", self.n_value, " (i.e ",self.n_value,"x_coord",self.n_value,")")
        print("TOTAL NUMBER OF RUNS: ", sum_runs)
        print()
        print("AVERAGE NUMBER OF RESTARTS: ", average_of_restart)
        print("AVERAGE NUMBER OF STEPS OF LAST RESTART ", avg_last_steps)
        print("AVERAGE NUMBER OF STEPS OF ALL RESTARTS", avg_total)
    
    
    # display_steepest_ascent_config -The report for steepest ascent hill climb  with sideways and without sideways move is displayed

    def display_steepest_ascent_config(self, result, first):
        
        sum_runs = result[0][0]
        succ = result[3][0]
        
        if succ:
            succ1=(succ/sum_runs)
            succ_rate = round((succ1)*100,2)
            success_steps = result[3][1]
            avearge_succ_steps = round(sum(success_steps)/succ, 2)
        else:
            succ_rate = success_steps = avearge_succ_steps = '-'
        
        fail = result[1][0]+result[2][0]
        
        if fail:
            fail1=(fail/sum_runs)
            fail_rate = round((fail1)*100,2)
            fail_steps = result[1][1]+result[2][1]
            fail_average_steps = round(sum(fail_steps)/fail,2)
        else:
            fail_rate = fail_steps = fail_average_steps = '-'
        
        num_flat_runs = result[1][0]
        
        print("\n\n"+first)
        border = ''
        l_f=len(first)
        for i in range(l_f): border+="="
        print(border)
        print("\n VALUE OF N: ", self.n_value, " (i.e ",self.n_value,"x",self.n_value,")")
        print("TOTAL NUMBER OF RUNS: ", sum_runs)
        print("\nSUCCESS RUNS: ", succ)
        print("SUCCESS RATES: ", succ_rate, "%")
        print("SUCCESS AVERAGE STEPS: ", avearge_succ_steps)
        print("\nNUMBER_OF_FAILURE RUNS: ", fail)
        print("FAILURE RATE: ", fail_rate, "%")
        print("FAILURE AVERAGE STEPS: ", fail_average_steps)
        print("\n\nFLAT LOCAL MINIMA: ", num_flat_runs)
        return


n_input = (int)(input("ENTER THE N VALUE GREATER THAN 3: "))

input_iterations = (int)(input("ENTER NUMBER OF ITERATIONS REQUIRED(MUST BE GREATER THAN OR EQUAL TO 1: "))
    
input_sideways = (int)(input("ENTER NUMBER OF MAXIMUM SIDEWAYS ALLOWED(MUST BE >=1): "))
   
if __name__ == "__main__":
    hill_climbing_analysis = Hill_Climbing_check(n_input, input_iterations, input_sideways)
    hill_climbing_analysis.explore()