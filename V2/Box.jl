#A set of 2D interacting particles. The interaction is modelled using a Morse Potential.
#The particles are initalized from auniform grid. 
using Plots
using JLD
plotlyjs()
##
function V_H(xl,yl)# Defing the interaction between particles to be harmonic 
    x = xl[1] - xl[2]
    y = yl[1] - yl[2]
    r = (x^2 + y^2)^0.5#Distance between the particles squared
    v = De*(1-exp(-a*(r-r_eq)))^2#Morse potential in 2 dimensions
    return v
end
##
function V_H_plot(xl,yl)# Defing the interaction between particles to be harmonic 
    x = xl
    y = yl
    r = (x^2 + y^2)^0.5#Distance between the particles squared
    v = De*(1-exp(-a*(r-r_eq)))^2#Morse potential in 2 dimensions
    return v
end
##
#Harmonic intearction is assumed between particles
function Fint_H(xl,yl)
    x = xl[1] - xl[2]
    y = yl[1] - yl[2]
    r = (x^2 + y^2)^0.5#Distance between the particles squared
    F = -2*De*(1-exp(-a*(r-r_eq)))*exp(-a*(r-r_eq))*a
    Fx, Fy = F*x/r, F*y/r
    return Fx,Fy
end
##
function Fx_H(xl,yl,k=K)
    x = xl
    y = yl
    r = (x^2 + y^2)^0.5#Distance between the particles squared
    F = -2*De*(1-exp(-a*(r-r_eq)))*exp(-a*(r-r_eq))*a
    Fx, Fy = F*x/r, F*y/r
    return Fx
end
function Fy_H(xl,yl,k=K)
    x = xl
    y = yl
    r = (x^2 + y^2)^0.5#Distance between the particles squared
    F = -2*De*(1-exp(-a*(r-r_eq)))*exp(-a*(r-r_eq))*a
    Fx, Fy = F*x/r, F*y/r
    return Fy
end
##
function Run(Xinit=0,Yinit=0,Vxinit=0,Vyinit=0;T=1:0.1:10,np=10)
    #INPUT: initial cordinates and velocities and the time axis and number of particles
    #RETURNS: The X,Y as a function of time. Total energy as a function of time
    lent = length(T)
    X = zeros(lent,np)
    Y = zeros(lent,np)
    E = zeros(lent)
    Temp = zeros(lent)#the temperature i.e, the average kinetic energy
    PE = zeros(lent)
    xnew = Xinit
    ynew = Yinit
    vxnew = Vxinit
    vynew = Vyinit
    for i in 1:lent
        #fxy = F.(xnew,[Xinit[3],Xinit[6]]ynew)
        #fxt = [ f[1] for f in fxy]
        #fyt = [ f[2] for f in fxy]
          #parameter for boundary
        xinit = copy(xnew)
        #
        fxt = zeros(np)
        fyt = zeros(np)
        for i in 1:length(xnew)-1
            for j in i+1:length(ynew)
                fxint,fyint = Fint_H([xnew[i],xnew[j]],[ynew[i],ynew[j]])
                fxt[i] += fxint
                fyt[i] += fyint
                fxt[j] += -fxint
                fyt[j] += -fyint
            end
        end
        xnew = xnew .+ vxnew*dt .+ 0.5*fxt*dt*dt./m
        ynew = ynew .+ vynew*dt .+ 0.5*fyt*dt*dt./m
        #fxyn = F.(xnew, ynew)
        #fxtn = [ f[1] for f in fxyn]
        #fytn = [ f[2] for f in fxyn]
        fxtn = zeros(np)
        fytn = zeros(np)
        PE_elec = 0 #Potential energy due to electronic interactions
        for i in 1:length(xnew)-1
            for j in i+1:length(ynew)
                fxint,fyint = Fint_H([xnew[i],xnew[j]],[ynew[i],ynew[j]])
                fxtn[i] += fxint
                fytn[i] += fyint
                fxtn[j] += -fxint
                fytn[j] += -fyint
                PE_elec += V_H([xnew[i],xnew[j]],[ynew[i],ynew[j]])
            end
        end
        vxnew = vxnew .+ 0.5*(fxt .+ fxtn)./m*dt
        vynew = vynew .+ 0.5*(fyt .+ fytn)./m*dt
        X[i,:] = xnew
        Y[i,:] = ynew
        PE[i] = PE_elec
        KE = sum( ( 0.5*m.*(vxnew.^2 + vynew.^2)) )
        E[i] =  KE + PE_elec
        Temp[i] = KE/np
        for i in 1:np
            #global box in x
            if (abs(xnew[i]) > bl)
                vxnew[i] = - vxnew[i]
            end
            if (abs(ynew[i]) > bl)
                vynew[i] = - vynew[i]
            end
        end
    end
    return X,Y,E,Temp,PE
end
##Interaction Parameters
np = 25
a = 0.5
De = 1
r_eq = 2.0
##Initialization Parameters
bl = 10 # half-box length
frames = 1000
##
print("Reading Initial Conditions")
file = jldopen("mydata.jld", "r")
    Xinit = read(file,"X") #.- 5 
    Yinit = read(file,"Y") #.- 5
close(file)
##
Yinit = Yinit .- (sum(Yinit)/length(Yinit))
Xinit = Xinit .- (sum(Xinit)/length(Xinit))
#print(Xinit)
np = length(Xinit)
dt =0.001
T = 0:dt:200
lt = length(T)
m = ones(np)
##
Vxinit = 5.0#0*(rand(np).-0.5)
Vyinit = 0.0#0.5#0*(rand(np).-0.5)
##
print("\nRunning Simulation") 
Xt,Yt,Et,Te,PE=Run(Xinit,Yinit,Vxinit,Vyinit;T=T,np=np)
print("\nEnergy Error:",maximum(Et)-minimum(Et))
##
print("\n Plotting")
plot(T,Et)
savefig("Total_Energy.png")
plot(T,Te)
savefig("Temperature.png")
##
#surface(-5:0.1:5,-5:0.1:5,V_H_plot)
#savefig("V.png")
##
N_step_frame =  Int(ceil(lt/frames))#No of steps per frame
print("\nAnimating\n")
plot(axis=([], false),ticks = false)
anim = @animate for j in 1:frames
    i = Int(ceil(j*lt/frames))
    plot([-bl,bl,bl,-bl,-bl],[-bl,-bl,bl,bl,-bl], color = :black,linewidth=3)
    if j > 10
        scatter!([Xt[i,:]],[Yt[i,:]],xlims=(-25,25),ylims=(-25,25),legend=false,color=:blue,markersize=2)
        plot!(Xt[(i-4*N_step_frame):i,:],Yt[(i-4*N_step_frame):i,:],legend=false,color=:blue,linewidth=1)
    else
        scatter!([Xt[i,:]],[Yt[i,:]],xlims=(-25,25),ylims=(-25,25),legend=false,color=:blue,markersize=2)
        plot!(Xt[1:i,:],Yt[1:i,:],legend=false,color=:blue,linewidth=1)
    end
end
gif(anim,"Morse_sim_box.gif", fps=30)
##
print("Done\n")