% initialize the staircase structure
sc = staircase_init(0.3,0.005);


nTrials = 300;

for t = 1:nTrials

acc = randsample([0 1 1 1],1)
[delta] = staircase_getDelta(sc);
sc = staircase_update(sc,delta,acc);

end

subplot(2,1,1);
plot(sc.acc)
subplot(2,1,2);
plot(sc.delta)
mean(sc.acc(100:nTrials))