%% loop 
[x, t] = readObj('../obj/horse.obj');
figure; subplot(221); trimesh(t, x(:,1), x(:,2), x(:,3), 'edgecolor', 'k'); axis equal; axis off; title('coarse mesh');

[x1, t1, T1] = loop(x, t);
subplot(222); trimesh(t1, x1(:,1), x1(:,2), x1(:,3), 'edgecolor', 'k'); axis equal; axis off; title('loop 1');

[x2, t2, T2] = loop(x1, t1);
subplot(223); trimesh(t2, x2(:,1), x2(:,2), x2(:,3), 'edgecolor', 'k'); axis equal; axis off; title('loop 2');

%% add noise
model_size = max(max(x)-min(x));
y2 = x2 + 0.006*model_size*rand(size(x2));
subplot(224); trimesh(t2, y2(:,1), y2(:,2), y2(:,3), 'edgecolor', 'k'); axis equal; axis off; title('loop 2 noise');

%% smooth 
s2 = smooth_loop(y2, {T2,T1}, 0.01*model_size); % 两层分解
s1 = smooth_loop(y2, {T2}, 0.01*model_size); % 一层分解
figure; subplot(131); trimesh(t2, y2(:,1), y2(:,2), y2(:,3), 'edgecolor', 'k'); axis equal; axis off; title('loop 2 noise');
subplot(132); trimesh(t2, s2(:,1), s2(:,2), s2(:,3), 'edgecolor', 'k'); axis equal; axis off; title('loop 2 smooth');
subplot(133); trimesh(t2, s1(:,1), s1(:,2), s1(:,3), 'edgecolor', 'k'); axis equal; axis off; title('loop 1 smooth');