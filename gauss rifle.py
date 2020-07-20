from vpython import *
k=10000000#彈性係數
b=0.003#空氣阻力係數(ok)
ball_radius=0.01195 #鐵球直徑(m)
ball_mass=0.0075 #鐵球質量(kg)
magnet_length=0.0098 #磁鐵(圓柱體)長(m)
magnet_mass=0.006 #磁鐵(圓柱體)質量(kg)
def AirResistance(obj):
    #print('air resistance:',-b*obj.v)
    return -b*obj.v
def MagneticForce(ball_idx,magnet_idx,balls):
    index_dist=abs(ball_idx-magnet_idx)-1
    dist=mag(balls[ball_idx].pos-balls[magnet_idx].real_pos)*1000#m->mm
    if index_dist ==0:
        f=14674288294.400 * (dist ** (-8.770))
    elif index_dist ==1:
        f=3560019737502910000000000000000000000000000*(dist**(-31.224))
    elif index_dist ==2:
        f=52745436146163100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000*dist**(-100.573)
    else:
        print('magnetic force is not set yet')

    return norm(balls[magnet_idx].real_pos-balls[ball_idx].pos)*f
def SpringForce(apos,aradius,bpos,bradius):
    if mag(apos-bpos)<aradius+bradius:
        return (aradius+bradius-mag(bpos-apos))*norm(apos-bpos)*k
    return vec(0,0,0)
class MagnetWithBalls(list):
    def __init__(self, position, radius, amount_of_frontBALL, amount_of_backBALL,*,first_ball=None,firstball_dist):
        self.radius=radius
        self.magnet=cylinder(radius=self.radius,color=color.red,index=amount_of_frontBALL+1,axis=vec(magnet_length,0,0),real_pos=vec(0,0,0),m=magnet_mass,v=vec(0,0,0),a=vec(0,0,0))

        ball=lambda: sphere(radius=self.radius,color=vec(0.5,0.5,0.5),m=ball_mass,v=vec(0,0,0),a=vec(0,0,0))
        if first_ball:
            self.isconnected_to_another=True
            self.start_calculate_idx=1
            self.first_ball=first_ball
        else:
            self.isconnected_to_another = False
            self.start_calculate_idx = 0
            self.first_ball=sphere(radius=self.radius,color=color.blue,v=vec(0.05,0,0),a=vec(0,0,0),m=ball_mass,opacity=0.1)
        super().__init__([self.first_ball]+[ball() for i in range(amount_of_frontBALL)]+[self.magnet]+[ball() for i in range(amount_of_backBALL)])

        self[-1].color=color.blue#最後一顆設成藍色

        for ball in self:
            if not hasattr(ball,'f_arrow'):
                ball.f_arrow = arrow(color=ball.color, shaftwidth=ball.radius / 2)

        self.setpos(position,firstball_dist,setFirstBall=not self.isconnected_to_another)
    def setpos(self,position,firstball_dist,setFirstBall=True):
        self.move_magnet(position)
        if setFirstBall:
            self.first_ball.pos = self.magnet.real_pos-vec(firstball_dist, 0, 0)
        for i in range(1,len(self)):
            try:
                self[i].pos=self.magnet.real_pos+vec(mag(self.magnet.axis)/2+(abs(i-self.magnet.index)*2-1)*self.radius,0,0)*((i-self.magnet.index)/abs(i-self.magnet.index))
            except ZeroDivisionError: #self[i]==self.magnet
                pass

    def move_magnet(self,pos):
        self.magnet.real_pos = pos
        self.magnet.pos = self.magnet.real_pos - vec(mag(self.magnet.axis) / 2, 0, 0)

    def cal_velocity(self):
        #歸零(但第一顆若是上一組最後一顆，則它的受力計算到一半，不能歸零!!!!!!)
        for i in range(self.start_calculate_idx,len(self)):
            self[i].f = vec(0, 0, 0)
        for i in range(len(self)):
            if i != self.magnet.index:
                magneticForceToTheBall = MagneticForce(i, self.magnet.index, self)
                self[i].f += magneticForceToTheBall
                self.magnet.f -= magneticForceToTheBall
            self[i].f += AirResistance(self[i])

            # leftball collision
            try:
                if i - 1 == self.magnet.index:
                    self[i].f += SpringForce(self[i].pos, self[i].radius, self[i - 1].real_pos,mag(self[i - 1].axis) / 2)
                elif i == self.magnet.index:
                    self[i].f += SpringForce(self[i].real_pos, mag(self[i].axis) / 2, self[i - 1].pos,self[i - 1].radius)
                else:
                    self[i].f += SpringForce(self[i].pos, self[i].radius, self[i - 1].pos, self[i - 1].radius)
            except IndexError:
                pass

            # rightball collision
            try:
                if i == self.magnet.index:
                    self[i].f += SpringForce(self[i].real_pos, mag(self[i].axis) / 2, self[i + 1].pos,self[i + 1].radius)
                elif i + 1 == self.magnet.index:
                    self[i].f += SpringForce(self[i].pos, self[i].radius, self[i + 1].real_pos,mag(self[i + 1].axis) / 2)
                else:
                    self[i].f += SpringForce(self[i].pos, self[i].radius, self[i + 1].pos, self[i + 1].radius)
            except IndexError:
                pass
    def update(self):
        # 若第一顆是上一組最後一顆，則它的位置已經計算過了!!!!!!)
        for i in range(self.start_calculate_idx,len(self)):
            self[i].a = self[i].f / self[i].m
            self[i].v += self[i].a * dt
            if i == self.magnet.index:
                self.move_magnet(self[i].real_pos + self[i].v * dt)
                self[i].f_arrow.pos = self[i].real_pos + vec(0, (i + 1) * self[i].radius * 2, 0)
            else:
                self[i].pos += self[i].v * dt
                self[i].f_arrow.pos = self[i].pos + vec(0, (i + 1) * self[i].radius * 2, 0)

            self[i].f_arrow.axis = self[i].f * 0.01


scene = canvas(width=1080, height=360,center = vec(0,0,0),background=color.white,range=0.1)

dt = 0.00001    #時間間隔
t = 0         #初始時間

gd= graph(width=1080, height=360,
      title='最末球速率',
      xtitle='t(s)', ytitle='v(m/s)')
balls=[]
for i in range(10):
    balls.append(MagnetWithBalls(vec(0.3*i,0,0),ball_radius/2,0,2,firstball_dist=0.05,first_ball=balls[i-1][-1] if i>0 else None))
    #balls[i].curve=gcurve(color=color.red)
last_balls=[balls[0][0]]+[i[-1] for i in balls]
for last_ball in last_balls:
    last_ball.curve=gcurve(color=color.red,interval=1000)
last_balls[0].curve.color=color.blue
max_v=gcurve(color=color.green,interval=1000)
text=label(color=color.black)

while True:
    rate(1/dt)
    t += dt            #累計時間
    text.text="time:%.6fs"%t

    for i in range(len(balls)):
        balls[i].cal_velocity()
    for i in range(len(balls)):
        balls[i].update()
        #balls[i].curve.plot(t,balls[i][-1].v.x)
    for i in last_balls:
        i.curve.plot(t,i.v.x)
    #for ball in set([ball for ball in ballset for ballset in balls]):
        #ball.f_arrow.pos=
    max_v.plot(t,max(*[ball.v.x for ball in last_balls]))
    scene.center=max([ball for ball in last_balls],key=lambda ball:ball.v.x).pos
    text.pos=scene.center+vec(0,ball_radius*2,0)