//////////////////////////
// Runge Kutta Fehlberg //
//////////////////////////

func rk45(from, h, to float64, y []float64, fn func(float64, []float64) []float64) [][]float64 {
  eps := 0.000001
	var parameters = len(y)

	// initialize 'outer slice'
	ySlice := [][]float64{}
	// initialize first 'inner slice'
	ySlice = append(ySlice,make([]float64, parameters+1))

	// fill with start values
	ySlice[0][0] = from
	for i := 0; i < parameters; i++ {
		ySlice[0][i+1] = y[i]
	}

	var k1 = make([]float64, parameters)
	var k2 = make([]float64, parameters)
	var k3 = make([]float64, parameters)
	var k4 = make([]float64, parameters)
	var k5 = make([]float64, parameters)
	var k6 = make([]float64, parameters)

	var yn = make([]float64, parameters)

	//
	var k2p = make([]float64, parameters)
	var k3p = make([]float64, parameters)
	var k4p = make([]float64, parameters)
	var k5p = make([]float64, parameters)
	var k6p = make([]float64, parameters)
	
  t := from;
	for  t < to {
	  step := len(ySlice)

		// initialize
		for value := 0; value < parameters; value++ {
			yn[value] = ySlice[step-1][value+1]
		}
    
    t1 := t
    
		// generate k1
		for value := 0; value < parameters; value++ {
			k1[value] = h * fn(t1, yn)[value]
		}
		

		// generate the parameter for k2
		for value := 0; value < parameters; value++ {
			k2p[value] = yn[value] + 0.25*k1[value]
		}
		
		t2 := t + 0.25 * h

		// generate k2
		for value := 0; value < parameters; value++ {
			k2[value] = h * fn(t2, k2p)[value]
		}

		// generate the parameter for k3
		for value := 0; value < parameters; value++ {
			k3p[value] = yn[value] + (3.0/32.0) * k1[value] + (9.0/32.0) * k2[value]
		}
		
		t3 := t + (3.0/8.0) * h

		// generate k3
		for value := 0; value < parameters; value++ {
			k3[value] = h * fn(t3, k3p)[value]
		}

		// generate the parameter for k4
		for value := 0; value < parameters; value++ {
			k4p[value] = yn[value] + (1932.0/2197.0) * k1[value] - (7200.0/2197) * k2[value] + (7296.0/2197.0) * k3[value]
		}

		t4 := t + (12.0/13.0) * h

		// generate k4
		for value := 0; value < parameters; value++ {
			k4[value] = h * fn(t4, k4p)[value]
		}
		
		// generate the parameter for k5
		for value := 0; value < parameters; value++ {
		  k5p[value] = yn[value] + (439.0/216.0) * k1[value] - 8.0 * k2[value] + (3680.0/513.0) * k3[value] - (845.0/4104.0) * k4[value]
		}

    t5 := t + h
    
    // generate k5
    for value := 0; value < parameters; value ++ {
      k5[value] = h * fn(t5,k5p)[value]
    }
    
    // generate the parameter for k6
    for value := 0; value < parameters; value ++ {
      k6p[value] = yn[value] - (8.0/27.0) * k1[value] + 2.0 * k2[value] - (3544.0/2565.0) * k3[value] + (1859.0/4104.0) * k4[value] - (11.0/40.0) * k5[value]
    }
    
    t6 := t + 0.5 * h
    
    // generate k6
    for value := 0; value < parameters; value ++ {
      k6[value] = h * fn(t6,k6p)[value]
    }
    
		
    // generate yk
    var yk = make([]float64, parameters)
    
    // generate the values for yk
    for value := 0; value < parameters; value ++ {
      yk[value] = yn[value] + (25.0/216.0) * k1[value] + (1408.0/2565.0) * k3[value] + (2197.0/4104.0) * k4[value] - (1.0/5.0) * k5[value]
    }
    
    
    //generate zk
    var zk = make([]float64, parameters)
    
    // generate the values for zk
    for value := 0; value < parameters; value ++ {
      zk[value] = yn[value] + (16.0/135.0) * k1[value] + (6656.0/12825.0) * k3[value] + (28561.0/56430.0) * k4[value] - (9.0/50.0) * k5[value] + (2.0/55.0) * k6[value]
    }
    
    
    // generate the residual
    L2_R := 0.0
    
    for value := 0; value < parameters; value ++ {
      L2_R += math.Pow(zk[value] - yk[value],2)
    }
    
    
    L2_R = math.Sqrt(L2_R)
    
    if L2_R <= eps {
      
      t = t + h
      
      // initialize 'inner slice'
		  ySlice = append(ySlice,make([]float64, parameters+1))
      
      ySlice[step][0] = t
      
      // Error small enough. Save yk slice to global
      for value := 0; value < parameters; value ++ {
        ySlice[step][value+1] = yk[value]
      }
      
    }
    
    // Find the scalar multiplier for the time step
    d := 0.84 * math.Pow(((eps*h)/L2_R),0.25)
    
    // Calculate the timestep
    h = d*h
    
    // If the timestep will put you over the final time, choose the difference as the timestep.
    if t + h > to {
      h = to-t
    }
    
    
    
	}
	return ySlice
}