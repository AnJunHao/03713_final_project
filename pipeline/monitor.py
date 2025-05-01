from datetime import datetime
import time
from pathlib import Path
def monitor_jobs(output_logs: list[Path],
                 error_logs: list[Path],
                 finish_signal: str = "Job finished") -> bool:
    """
    Monitor job status by checking output and error logs.
    
    Args:
        output_logs: List of paths to output log files
        error_logs: List of paths to error log files
        
    Returns:
        True when all jobs have completed successfully
    """
    job_statuses = {str(log): "QUEUED" for log in output_logs}
    all_complete = False
    start_time = datetime.now()
    
    # Print header once
    print("Monitoring job status...")
    
    # Initialize the status lines but don't print them yet
    status_lines = []
    for _ in range(len(output_logs)):
        status_lines.append("")
    
    # Add an extra line for elapsed time
    status_lines.append("")
    
    while not all_complete:
        all_complete = True
        
        # Clear previous status lines if any exist
        if status_lines[0]:  # If we've already printed status lines
            for _ in range(len(output_logs) + 1):  # +1 for elapsed time
                print("\033[A\033[K", end="")
        
        # Check and update each job's status
        for i, (out_log, err_log) in enumerate(zip(output_logs, error_logs)):
            job_name = out_log.name.replace(".out.txt", "")
            status_icon = "📋"
            status_text = "QUEUED"
            detail_text = ""
            
            # Check if logs exist
            if not out_log.exists() and not err_log.exists():
                job_statuses[str(out_log)] = "QUEUED"
                all_complete = False
            
            # Check if error log has content
            elif err_log.exists() and err_log.stat().st_size > 0:
                with open(err_log, 'r') as f:
                    err_content = f.read().strip()
                    if err_content:
                        job_statuses[str(out_log)] = "ERROR"
                        status_icon = "❌"
                        status_text = "ERROR"
                        # Get first line or first 50 chars of error
                        err_lines = err_content.splitlines()
                        detail_text = err_lines[-1][:50] + ("..." if len(err_lines[-1]) > 50 else "")
            
            # Check if output log exists and job is complete
            elif out_log.exists():
                if job_statuses[str(out_log)] == "QUEUED":
                    job_statuses[str(out_log)] = "RUNNING"
                
                # Check if job is complete by reading last line
                with open(out_log, 'r') as f:
                    content = f.read().strip()
                    if content:
                        lines = content.splitlines()
                        last_line = lines[-1] if lines else ""
                        
                        # Get last line of output for detail text
                        if lines:
                            detail_text = lines[-1][:50] + ("..." if len(lines[-1]) > 50 else "")
                        
                        if last_line == finish_signal:
                            job_statuses[str(out_log)] = "COMPLETED"
                            status_icon = "✅"
                            status_text = "COMPLETED"
                        else:
                            job_statuses[str(out_log)] = "RUNNING"
                            status_icon = "🏃"
                            status_text = "RUNNING"
                            all_complete = False
                    else:
                        job_statuses[str(out_log)] = "RUNNING"
                        status_icon = "🏃"
                        status_text = "RUNNING"
                        all_complete = False
            else:
                all_complete = False
            
            # Format and store the status line
            status_lines[i] = f"{status_icon} {job_name} {status_text}: {detail_text}"
            
            # Print current status
            print(status_lines[i])
        
        # Calculate and display elapsed time
        elapsed = datetime.now() - start_time
        hours, remainder = divmod(elapsed.total_seconds(), 3600)
        minutes, seconds = divmod(remainder, 60)
        elapsed_str = f"⏱️ Elapsed time: {int(hours):02d}:{int(minutes):02d}:{int(seconds):02d}"
        status_lines[-1] = elapsed_str
        print(elapsed_str)
        
        if not all_complete:
            time.sleep(1)  # Refresh every second
    
    # Count completed and error jobs
    completed_jobs = sum(1 for status in job_statuses.values() if status == "COMPLETED")
    error_jobs = sum(1 for status in job_statuses.values() if status == "ERROR")
    total_jobs = len(job_statuses)
    
    # Print final status with appropriate emoji
    if error_jobs == 0:
        print(f"🎉 All {total_jobs} jobs completed successfully in {elapsed_str[17:]}")
    elif completed_jobs == 0:
        print(f"❌ All {total_jobs} jobs failed")
        # Print error log locations for all failed jobs
        print("\nError logs:")
        for err_log in error_logs:
            if err_log.exists() and err_log.stat().st_size > 0:
                print(f"  - {err_log.absolute()}")
                # Show first 3 lines of error log
                with open(err_log, 'r') as f:
                    first_lines = [line.strip() for line in f.readlines()[:3]]
                    for line in first_lines:
                        print(f"      {line}")
    else:
        print(f"⚠️ {completed_jobs}/{total_jobs} jobs completed, {error_jobs} jobs failed")
        # Print error log locations for failed jobs
        print("\nError logs:")
        for i, (out_log, err_log) in enumerate(zip(output_logs, error_logs)):
            if job_statuses[str(out_log)] == "ERROR" and err_log.exists() and err_log.stat().st_size > 0:
                print(f"  - {err_log.absolute()}")
                # Show first 3 lines of error log
                with open(err_log, 'r') as f:
                    first_lines = [line.strip() for line in f.readlines()[:3]]
                    for line in first_lines:
                        print(f"      {line}")
    
    return error_jobs == 0