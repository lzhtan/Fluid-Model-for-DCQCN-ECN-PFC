
%{
(1) Program: Simulate packet loss of the bottleneck switch in the following topology
(2) Topology：
            Server1 <-> S0(S0) <-> S1 <-> S2 <-> S1 <-> S0(ToR) <-> Server2
            where，S2 only enables ECN, not PFC; both S0 and S1 enable ECN and PFC.
%}
function sol = dcqcn2()
    clc;
    clear all;
    close all;
    
%%%%%%% 1. Network parameters %%%%%%%%

    global packetSize; %Packet size, DCQCN paper is set to 1KB (8000 bits)
    packetSize = 8e3; % in bits!
    
    global C; %Bottleneck switch capacity, DCQCN paper is set to 40 Gbps, corresponding to 5 million packets
    C = 40e9/packetSize;
    
    global numFlows;  % Number of streams
    
    global initVal; %Store the state set of each flow, including the current rate, target rate, and current alpha value
    
    global lossNumPakets; %丢包数量
    lossNumPackets = 0;

%%%%%%% 2. DCQCN parameters %%%%%%%%

    global Rai; %In the fastrecovery stage, a constant (step size) is added to the target rate, and the DCQCN paper is set to 40 Mbps for Rt = Rt + Rai
    Rai = 40e6/packetSize;

    global B; %Bytecounter, DCQCN byte counter, DCQCN paper is set to 10 MB, corresponding to 10000 packets
    B = 10e6*8/packetSize;
    
    global timer; %Rate increase timer, DCQCN paper set to 55us
    timer = 55e-6;
    
    global F; %Number of fastrecovery counters, the DCQCN paper is set to 5
    F = 5;
    
    global g; %Alpha adjustment coefficient, DCQCN papers are set to 1/16 and 1/256
    g = 1/256;
    
    global tau; %Response cycle of the receiving end, DCQCN paper is set to 50 us, and the receiving end can generate at most 1 CNP within 50us
    tau = 50e-6;
    
    global taustar; %Feedback delay, the paper is set to 50us, and the simulation is set to 85 us, which can be understood as a delay parameter related to rtt and the response interval of the receiving end
    
    global tauprime; %If the sending end does not receive cnp within 55s, then reduce the alpha value to alpha = (1-g)alpha
    tauprime = 55e-6;
    
%%%%%%% 3. ECN/PFC parameters %%%%%%%%

    global queueLength; %Queue length, set to 1000 K, corresponding to 1000 packets
    queueLength = 10e5*8/ packetSize;
    
    global Kmin; %ECN low watermark, set to 5 KB, corresponding to 5 packets
    Kmin = 5e3*8/ packetSize;
    
    global Kmax; %ECN high watermark, set to 200 KB, corresponding to 200 packets
    Kmax = 200e3*8/ packetSize;
    
    global pmax; %ECN maximum marking probability, set to 10%
    pmax = 5e-1;
    
    global pfcThreshold; %PFC trigger threshold, set to 300 K, corresponding to 300 packets
    pfcThreshold = 800e3*8/ packetSize;

    global pfcPauseTime; %Pause duration after PFC is triggered. The default configuration of the switch is usually to forward 65535*512bit packets, which is about 838 us in a 40G network environment
    pfcPauseTime = 65535*512/40e9;
    
    global pfcStartTime; %PFC trigger time starting point, last pfcPauseTime
    pfcStartTime = 0;
    
   
    
%%%%%%% 4. Parameters for simulation %%%%%%%%
    global numCalls; %Simulation rounds
    numCalls = 0;
 
    sim_step = 10e-6; %Simulation time step, set to 10us
    
    options = ddeset('MaxStep', sim_step, 'RelTol', 1e-2, 'AbsTol', 1e-4);
    %Create or change the delay differential equation options structure, MaxStep is the maximum step size of 50 ms, RelTol is the relative error tolerance of 1%, and AbsTol is the absolute error tolerance of 1%.    
    sim_length = 50e-3;%The total time length of the simulation, set to 50ms
    
    % !!!!!!
    % All quantities (rates, buffers, queues etc.) are specified in units of packet size. The reason
    % is that the probability calculations are per packet.b
    % !!!
   

    for taustar = [85e-6]  % vary the feedback delay
        for numFlows = [8]
            % Initial conditions: (single column matrix)
            %
            % 1: rc1
            % 2: rt1
            % 3: alpha1
            % 4: rc2
            % 5: rt2
            % 6: alpha2
            % 7: ...
            % 3*numFlows+1: queue2 -> S1
            % 3*numFlows+2: queue1 -> S0
            % 3*numFlows+3: packet loss

            % set initial and traget rate of all flows to C/N, and alpha to 1
            %Specify three state parameters for each flow: current rate/target rate/alpha value, which is an array of 3*numFlows rows and 1 column, and the 3rd*numFlows+1 row is the queue length, initialized to 0            
            initVal = zeros(3*numFlows + 3, 1); 
            initVal(1:3:3*numFlows,1) = C; %Initialize the current rate
            initVal(2:3:3*numFlows,1) = C; %Initialize target rate
            initVal(3:3:3*numFlows,1) = 1; %Initialize the alpha value
        
            % solve.
            sol = dde23(@DCQCNModel, taustar, initVal, [0, sim_length], options);
            %This is the most time consuming part of the whole program
            %Differential equation solver dde23, DCQCNModel is a differential equation; taustar is a time delay of 85us; initVal is an initial value; options are previously defined parameters, including the total simulation time range and tolerable error.            
            % Extract solution and write to file.
            t = sol.x; 
            queueS1 = sol.y(end-2,:); %S1 queue length changes
            queueS0 = sol.y(end-1,:); %S0 queue length changes
            rates = sol.y(1:3:end-3,:); %Actual rate per flow
            lossPacketNum = sol.y(end,:); %Changes in packet loss
            [utilization, err] = Utilization(t, rates, queueS0, C);  
            fprintf('utilization: flows=%d util=%f err=%d\n', numFlows, utilization, err);              

            prefix = 'unstable';
            fileName = sprintf('%s.%d.%d.dat', prefix, numFlows, taustar*1e6);
            
            % when writing to file, rate is in Gbps, and queue is in KB.
            %dlmwrite(fileName,[t',rates'.*packetSize/1e9, queueS0'.*packetSize/8e3], '\t');
            fprintf ('S1 median queue=%f\n', median(queueS1));
            fprintf ('S0 median queue=%f\n', median(queueS0));

            PlotSol(t, queueS0, queueS1,lossPacketNum, rates, sim_length);
            
        end
        fclose('all');
    end
end
    
function dx = DCQCNModel(t,x,lag_matrix) %DCQCN model, t is the simulation time, x is the initial state, lag_matrix stores the current state
    global Kmax;
    global Kmin;
    global numFlows;
    global numCalls;
    global lossNumPakcets;
    
    % matrix x:
    % 1: rc1
    % 2: rt1
    % 3: alpha1
    % 4: rc2
    % 5: rt2
    % 6: alpha2
    % ...
    % 3*numFlows+1: queue2 -> S1
    % 3*numFlows+2: queue1 -> S0
    
    % lag matrix: 
    % (:,1) is t-t' for flow 1
    % (:,2) is t-t' for flow 2
    % ....

    dx = zeros(3*numFlows+3,1);%Store the current rate, target rate, and alpha value change of the current round
    
    %
    % marking probability
    %
    
    p = CalculateP(t,lag_matrix(end-1,1),Kmin,Kmax);%Calculation of marking probabilities
    %
    
    %
    % rates and alpha
    %
    for i = 1:3:3*numFlows 
        % The model cannot correctly handle the case of prevRC = 0. So we stay at 0. 
        if (x(i) == 0 && lag_matrix(i, 1) == 0) %If the current rate of each flow is 0 and the rate in the lag is also 0, set the current rate and target rate change to 0
            dx(i) = 0; 
            dx(i+1) = 0; 
        else
            [a,b,c,d,e] = IntermediateTerms(p, lag_matrix(i, 1), x(i), t, i);
            % rc
            dx(i) =  RCDelta(x(i), x(i+1), x(i+2), lag_matrix(i, 1), a, b, d);
            % rt
            dx(i+1) = RTDelta(x(i), x(i+1), lag_matrix(i, 1), a, c, e);
        end
        % alpha
        dx(i+2) = AlphaDelta(x(i+2), lag_matrix(i, 1), p);
    end
    
    %
    % Queue
    %
    rates = x(1:3:3*numFlows);
    [dx(end-2), dx(end-1), dx(end)] = QueueDelta(x(end-2), x(end-1), rates, t);

    
     numCalls = numCalls +1;
     if (mod(numCalls, 10000) == 0)
         fprintf ('%g %d %f\n', t, numCalls, p);
     end
    
end

%prevRC is Rc(t-tau*) in the Fluid model
%Current rate change RcDelta, the input is the current actual rate, current target rate, current Alpha value, Rc(t-tau*), temporary values in the three groups of Fluid models
function rcDelta = RCDelta(currRC, currRT, currAlpha, prevRC, a, b, d)
    global tau;
    global C;
   
    rcDelta = -1*currRC*currAlpha*a/(2*tau) + (currRT-currRC)*prevRC*b/2 + (currRT-currRC)*prevRC*d/2;
    
    %Rate cannot exceed C.
    if (currRC >= C && rcDelta > 0)
        rcDelta = 0;
    end
end

%Target rate change RtDelta, the input is the current actual rate, current target rate, Rc(t-tau*), temporary values in the three groups of Fluid models
function rtDelta = RTDelta(currRC, currRT, prevRC, a, c, e)
    global tau;
    global Rai;
    global C;
    
    rtDelta = -1*(currRT-currRC)*a/(tau) + Rai*prevRC*c + Rai*prevRC*e;
    
    % rate cannot exceed C.
    if (currRT >= C && rtDelta > 0)
        rtDelta = 0;
    end 
end


function alphaDelta = AlphaDelta(currentAlpha, prevRC, p)
    global tauprime;
    global g;
    alphaDelta = g * ((1-(1-p)^(tauprime*prevRC))-currentAlpha) / tauprime;
end



function [S1queueDelta, S0queueDelta, lossNumPakcetsDelta] = QueueDelta(S1currentQueue, S0currentQueue, rates, t)
    global C;
    global pfcStartTime;
    global pfcPauseTime;
    global queueLength;
    global pfcThreshold;
    % Judging that if PFC is triggered, S1 can only enter but not enter, and S0 can only exit but not enter
    if pfcStartTime ~= 0
        if (t - pfcStartTime) < pfcPauseTime
            S1queueDelta = sum(rates); %After PFC is triggered, the S1 switch queue continues to increase
            S0queueDelta = - C;
        elseif (t - pfcStartTime) >= pfcPauseTime
            if S0currentQueue >= pfcThreshold
                pfcStartTime = t;
                S1queueDelta = sum(rates); %After PFC is triggered, the S1 switch queue continues to increase
                S0queueDelta = - C;
            elseif S0currentQueue < pfcThreshold
                S0queueDelta = sum(rates) - C;
                S1queueDelta = sum(rates) - C;
            end
        end
    end
    
    if pfcStartTime == 0
        if S0currentQueue < pfcThreshold
            S0queueDelta = sum(rates) - C;
            S1queueDelta = sum(rates) - C;
        elseif S0currentQueue >= pfcThreshold
            pfcStartTime = t;
            S1queueDelta = sum(rates); %After PFC is triggered, the S1 switch queue continues to increase
            S0queueDelta = - C;
        end
    end
    
    if (S0currentQueue <= 0) && (S0queueDelta < 0)
        S0queueDelta = max(0, S0queueDelta);
    end
    
    if (S1currentQueue <= 0) && (S1queueDelta < 0)
        S1queueDelta = max(0, S1queueDelta);
    end
    
    if ((S1currentQueue >= queueLength) && (S1queueDelta > 0)) % If the queue length exceeds 1000 packets, the queue will no longer increase
        S1queueDelta = 0;
        lossNumPakcetsDelta = sum(rates); %The number of packets dropped
    else
        lossNumPakcetsDelta = 0;
    end
   
    if (S0currentQueue >= queueLength) 
        S0queueDelta = 0;
    end
end

function p = CalculateP(t, q, kmin, kmax)
    global pmax;
    if (t >= 0)
        if q <= kmin
            p = 0;
        else if q <= kmax
                p = (q-kmin)/(kmax-kmin)*pmax;
            else % q > lmax
                p = 1;
            end
        end
    else 
        p = 1;
    end 
end

% Input：
% mark probability, prevRC, current rate currRC

% Output:
% a corresponds to (1-(1-p(t-tau*))^{tau*RcDleta}) in the Fluid model
% b corresponds to p(t-tau*)/[(1-p(t-tau*))^{-B}-1] in the Fluid model, used for RcDelta
% c corresponds to (1-p(t-tau*))^{FB}*p(t-tau*)/[(1-p(t-tau*))^{-B}-1] in the Fluid model , for RtDelta
% d corresponds to the last fraction of RtDelta in the Fluid model
% e corresponds to the last fraction of RcDelta in the Fluid model

function [a, b, c, d, e] = IntermediateTerms(p, prevRC, currRC, t, i)
    
    global tau;
    global B;
    global F;
    global timer;   
    if p == 0 
            a = 0;
            b = 1/B;   
            c = b;
            if (prevRC == 0)
                d = 0;
            else 
                d = 1/(timer*prevRC);
            end
            e = d;
        else if p == 1 
                a = 1;
                b = 0;
                c = 0;
                d = 0;
                e = 0;
            else    
                a = 1-(1-p)^(tau*prevRC);   
                b = p/((1-p)^(-B)-1);
                c = b*((1-p)^(F*B));
                d = p/((1-p)^(-timer*prevRC)-1);
                if (isinf(d))
                    d  = 1/(timer*prevRC);
                    e = d;
                else
                    e = d*((1-p)^(F*timer*prevRC));
                end
            end
    end
    if (isnan(a) || isnan(b) || isnan(c) || isnan(d) || isinf(d) || isinf(e))
        fprintf (' *************** ERROR ***********************\n');
        fprintf ('%i ', [p prevRC currRC t i a b c d e]);
        fprintf ('\n');
        a = 0; 
        b = 1/B;
        if (prevRC == 0)
           d = 0;
        else 
           d = 1/(timer*prevRC);
        end
        e = d;
    end
end

function [u, err] = Utilization (t, rates, q, C)
    sent = 0;
    tmin = t(1,1);
    tmax = t(1,end);
    max = C * (tmax - tmin);
    err = 0;
    for tindex = 1:(size(t, 2)-1)
        ratesum = 0;
        if (q(tindex) > 1)
            ratesum = C;
        else 
            for flow = 1:size(rates, 1)
                ratesum = ratesum + rates(flow, tindex);
            end
            if (ratesum > C)
                ratesum = C;
                err = err + 1;
            end
        end
        sent = sent + ratesum * (t(1, tindex+1) - t(1, tindex) );
    end
    u = sent/max;
end

function PlotSol(t, queueS0, queueS1, lossPacketNum, rates, sim_length)
    global packetSize;
    figure
    subplot(4,1,1);
    Figure1=rates'*packetSize/1e9
    plot(t,rates'*packetSize/1e9,'LineWidth',2);
    hold on
    axis([0 sim_length 0 1.1*max(max(rates'*packetSize/1e9))])
    xlabel('Time (seconds)')
    ylabel('Throughput (Gbps)')
    title('Changes in per-flow throughput');

    
    subplot(4,1,2);
    plot(t,queueS0,'LineWidth',2)
    Figure2=queueS0';
    hold on
    axis([0 sim_length 0 1.1*max(queueS0)])
    xlabel('Time (seconds)')
    ylabel('Queue (Packets)')
    title('Queue changes of switch S0（with ECN&PFC）');

    
    subplot(4,1,3);
    plot(t,queueS1,'LineWidth',2)
    hold on
    Figure3=queueS1';
    axis([0 sim_length 0 1.1*max(queueS1)])
    xlabel('Time (seconds)')
    ylabel('Queue (Packets)')
    title('Queue changes of switch S1（with ECN&PFC）');
    
    subplot(4,1,4);
    plot(t,lossPacketNum,'LineWidth',2)
    hold on
    Figure4=lossPacketNum'
    axis([0 sim_length 0 max(1, 1.1*max(lossPacketNum))])
    xlabel('Time (seconds)')
    ylabel('Queue (Packets)')
    title(' Packet loss changes of switch S1');

    

end
