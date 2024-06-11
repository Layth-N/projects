int trigPin = 2;    // Trigger
int echoPin1 = 8;    // Echo
long duration1, cm1;
int echoPin2 = 9;    // Echo
long duration2, cm2;
int echoPin3 = 10;    // Echo
long duration3, cm3;
int echoPin4 = 11;    // Echo
long duration4, cm4;

//BackWheel Drive

//Channel pins
#define CH1 3
#define CH2 5
#define CH9 6

//Relays Pins
int F1 = A0;
int F2 = A1;
int BA1 = A2;
int BA2 = A3;
int S1 = A4;
int S2 = A5;
int UV1 = 4;
int UV2 = 7;


// Read the number of a given channel and convert to the range provided.
// If the channel is off, return the default value
int readChannel(int channelInput, int minLimit, int maxLimit, int defaultValue){
  int ch = pulseIn(channelInput, HIGH, 30000);
  if (ch < 100) return defaultValue;
  return map(ch, 1000, 2000, minLimit, maxLimit);
}


void setup(){

  pinMode(trigPin, OUTPUT);
  pinMode(echoPin1, INPUT);
  pinMode(echoPin2, INPUT);
  pinMode(echoPin3, INPUT);
  pinMode(echoPin4, INPUT);
  
  pinMode(CH1, INPUT);
  pinMode(CH2, INPUT);
  pinMode(CH9, INPUT);
  
  pinMode(F1, OUTPUT);
  pinMode(F2, OUTPUT);
  pinMode(BA1, OUTPUT);
  pinMode(BA2, OUTPUT);
  pinMode(S1, OUTPUT);
  pinMode(S2, OUTPUT);

  //default relays state
  digitalWrite (F1,HIGH);
  digitalWrite (F2,HIGH);
  digitalWrite (BA1,HIGH);
  digitalWrite (BA2,HIGH);
  digitalWrite (S1,HIGH);
  digitalWrite (S2,HIGH);
  digitalWrite (UV1,HIGH);
  digitalWrite (UV2,HIGH);
}

int ch1Value, ch2Value, ch9Value;

void loop() {

 digitalWrite(trigPin, LOW);
 delayMicroseconds(5);
 digitalWrite(trigPin, HIGH);
 delayMicroseconds(10);
 digitalWrite(trigPin, LOW);
 
 // Read the signal from the sensor: a HIGH pulse whose
 // duration is the time (in microseconds) from the sending
 // of the ping to the reception of its echo off of an object.
 duration1 = pulseIn(echoPin1, HIGH);
 duration2 = pulseIn(echoPin2, HIGH);
 duration3 = pulseIn(echoPin3, HIGH);
 duration4 = pulseIn(echoPin4, HIGH);
 
 // Convert the time into a distance
 cm1 = (duration1/2) / 29.1;     // Divide by 29.1 or multiply by 0.0343
 cm2 = (duration2/2) / 29.1;
 cm3 = (duration3/2) / 29.1;
 cm4 = (duration4/2) / 29.1;

if (cm1 < 5 && cm2 < 5 && cm3 > 5 && cm4 > 5){

  ch1Value = readChannel(CH1, -100, 100, 0);
  ch2Value = readChannel(CH2, -100, 100, 0);
  ch9Value = readChannel(CH9, -100, 100, 0);

//Going Backward at LOW Speed
  digitalWrite (F1,HIGH);
  digitalWrite (F2,HIGH);
  digitalWrite (BA1,LOW);
  digitalWrite (BA2,LOW);
  digitalWrite (S1,HIGH);
  digitalWrite (S2,HIGH);
}

if (cm1 > 5 && cm2 > 5 && cm3 < 5 && cm4 < 5){

  ch1Value = readChannel(CH1, -100, 100, 0);
  ch2Value = readChannel(CH2, -100, 100, 0);
  ch9Value = readChannel(CH9, -100, 100, 0);
//Going Forward at LOW Speed
  digitalWrite (F1,LOW);
  digitalWrite (F2,LOW);
  digitalWrite (BA1,HIGH);
  digitalWrite (BA2,HIGH);
  digitalWrite (S1,HIGH);
  digitalWrite (S2,HIGH);
}

if (cm1 > 5 && cm2 > 5 && cm3 > 5 && cm4 > 5){

  ch1Value = readChannel(CH1, -100, 100, 0);
  ch2Value = readChannel(CH2, -100, 100, 0);
  ch9Value = readChannel(CH9, -100, 100, 0);

//Going Forward at HIGH Speed
if (ch2Value>60 && ch1Value<10 && ch1Value>-10) {
  digitalWrite (F1,LOW);
  digitalWrite (F2,LOW);
  digitalWrite (BA1,HIGH);
  digitalWrite (BA2,HIGH);
  digitalWrite (S1,LOW);
  digitalWrite (S2,LOW);
}

//Going Forward at LOW Speed
if (ch2Value<60 && ch2Value>10 && ch1Value<10 && ch1Value>-10) {
  digitalWrite (F1,LOW);
  digitalWrite (F2,LOW);
  digitalWrite (BA1,HIGH);
  digitalWrite (BA2,HIGH);
  digitalWrite (S1,HIGH);
  digitalWrite (S2,HIGH);
}

//Going Backward at HIGH Speed
if (ch2Value<-70 && ch1Value<10 && ch1Value>-10) {
  digitalWrite (F1,HIGH);
  digitalWrite (F2,HIGH);
  digitalWrite (BA1,LOW);
  digitalWrite (BA2,LOW);
  digitalWrite (S1,LOW);
  digitalWrite (S2,LOW);
}

//Going Backward at LOW Speed
if (ch2Value>-70 && ch2Value<-10 && ch1Value<10 && ch1Value>-10) {
  digitalWrite (F1,HIGH);
  digitalWrite (F2,HIGH);
  digitalWrite (BA1,LOW);
  digitalWrite (BA2,LOW);
  digitalWrite (S1,HIGH);
  digitalWrite (S2,HIGH);
}

//Turnin Right Forward At HIGH Speed
if (ch2Value>60 && ch1Value>50) {
  digitalWrite (F1,LOW);
  digitalWrite (F2,LOW);
  digitalWrite (BA1,HIGH);
  digitalWrite (BA2,HIGH);
  digitalWrite (S1,LOW);
  digitalWrite (S2,HIGH);
}

//Turnin Right Forward At LOW Speed
if (ch2Value<60 && ch2Value>10 && ch1Value>50) {
  digitalWrite (F1,LOW);
  digitalWrite (F2,HIGH);
  digitalWrite (BA1,HIGH);
  digitalWrite (BA2,HIGH);
  digitalWrite (S1,HIGH);
  digitalWrite (S2,HIGH);
}

//Turnin Right Backward At HIGH Speed
if (ch2Value<-70 && ch1Value>50) {
  digitalWrite (F1,HIGH);
  digitalWrite (F2,HIGH);
  digitalWrite (BA1,LOW);
  digitalWrite (BA2,LOW);
  digitalWrite (S1,LOW);
  digitalWrite (S2,HIGH);
}

//Turnin Right Backward At LOW Speed
if (ch2Value>-70 && ch2Value<-10 && ch1Value>50) {
  digitalWrite (F1,HIGH);
  digitalWrite (F2,HIGH);
  digitalWrite (BA1,LOW);
  digitalWrite (BA2,HIGH);
  digitalWrite (S1,HIGH);
  digitalWrite (S2,HIGH);
}

//Turnin Left Forward At HIGH Speed
if (ch2Value>60 && ch1Value<-50) {
  digitalWrite (F1,LOW);
  digitalWrite (F2,LOW);
  digitalWrite (BA1,HIGH);
  digitalWrite (BA2,HIGH);
  digitalWrite (S1,HIGH);
  digitalWrite (S2,LOW);
}

//Turnin Left Forward At LOW Speed
if (ch2Value<60 && ch2Value>10 && ch1Value<-50) {
  digitalWrite (F1,HIGH);
  digitalWrite (F2,LOW);
  digitalWrite (BA1,HIGH);
  digitalWrite (BA2,HIGH);
  digitalWrite (S1,HIGH);
  digitalWrite (S2,HIGH);
}

//Turnin Left Backward At HIGH Speed
if (ch2Value<-70 && ch1Value<-50) {
  digitalWrite (F1,HIGH);
  digitalWrite (F2,HIGH);
  digitalWrite (BA1,LOW);
  digitalWrite (BA2,LOW);
  digitalWrite (S1,HIGH);
  digitalWrite (S2,LOW);
}

//Turnin Left Backward At LOW Speed
if (ch2Value>-70 && ch2Value<-10 && ch1Value<-50) {
  digitalWrite (F1,HIGH);
  digitalWrite (F2,HIGH);
  digitalWrite (BA1,HIGH);
  digitalWrite (BA2,LOW);
  digitalWrite (S1,HIGH);
  digitalWrite (S2,HIGH);
}

//Turning Right Inplace at HIGH Speed
if (ch1Value>50 && ch2Value<10 && ch2Value>-10) {
  digitalWrite (F1,LOW);
  digitalWrite (F2,HIGH);
  digitalWrite (BA1,HIGH);
  digitalWrite (BA2,LOW);
  digitalWrite (S1,LOW);
  digitalWrite (S2,LOW);
}

//Turning Right Inplace at Low Speed***
if (ch1Value<50 && ch1Value>10 && ch2Value<10 && ch2Value>-10) {
  digitalWrite (F1,LOW);
  digitalWrite (F2,HIGH);
  digitalWrite (BA1,HIGH);
  digitalWrite (BA2,LOW);
  digitalWrite (S1,HIGH);
  digitalWrite (S2,HIGH);
}

//Turning Left Inplace at HIGH Speed
if (ch1Value<-50 && ch2Value<10 && ch2Value>-10) {
  digitalWrite (F1,HIGH);
  digitalWrite (F2,LOW);
  digitalWrite (BA1,LOW);
  digitalWrite (BA2,HIGH);
  digitalWrite (S1,LOW);
  digitalWrite (S2,LOW);
}

//Turning Left Inplace at LOW Speed
if (ch1Value>-50 && ch1Value<-10 && ch2Value<10 && ch2Value>-10) {
  digitalWrite (F1,HIGH);
  digitalWrite (F2,LOW);
  digitalWrite (BA1,LOW);
  digitalWrite (BA2,HIGH);
  digitalWrite (S1,HIGH);
  digitalWrite (S2,HIGH);
}

//Back to Default State
if (ch1Value<10 && ch1Value>-10 && ch2Value<10 && ch2Value>-10) {
  digitalWrite (F1,HIGH);
  digitalWrite (F2,HIGH);
  digitalWrite (BA1,HIGH);
  digitalWrite (BA2,HIGH);
  digitalWrite (S1,HIGH);
  digitalWrite (S2,HIGH);
}
}

//ALL UV OFF
if (ch9Value<-50) {
  digitalWrite (UV1,HIGH);
  digitalWrite (UV2,HIGH);
}

//ALL UV OFF
//if (ch9Value==-101) {
 // digitalWrite (UV1,HIGH);
  //digitalWrite (UV2,HIGH);
//}

//5 UV ON
if (ch9Value<50 && ch9Value>-50) {
  digitalWrite (UV1,LOW);
  digitalWrite (UV2,HIGH);
}

//5 UV ON
//if (ch9Value==-3) {
 // digitalWrite (UV1,LOW);
 // digitalWrite (UV2,HIGH);
//}

//ALL UV ON
if (ch9Value>50) {
  digitalWrite (UV1,LOW);
  digitalWrite (UV2,LOW);
}

//ALL UV ON
//if (ch9Value==97) {
 // digitalWrite (UV1,LOW);
 // digitalWrite (UV2,LOW);
//}

delay(250);

}
