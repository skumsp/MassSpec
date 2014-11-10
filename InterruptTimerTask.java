package massspec;

import java.util.TimerTask;


/*
 * A TimerTask that interrupts the specified thread when run.
 */
public class InterruptTimerTask extends TimerTask {

    private Thread theTread;

    public InterruptTimerTask(Thread theTread) {
        this.theTread = theTread;
    }

    @Override
    public void run() {
        theTread.interrupt();
    }

}