custom_css = """
.main-title {
    color: #2c3e50;
    margin-bottom: 2rem;
    text-align: center;
    padding: 1rem;
    background: linear-gradient(to right, #f8f9fa, #e9ecef);
    border-radius: 8px;
}

.upload-section {
    background-color: #ffffff;
    border-radius: 10px;
    padding: 2rem;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
}

.analysis-card {
    background-color: #ffffff;
    border-radius: 10px;
    padding: 1.5rem;
    margin: 1rem 0;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
}

.consent-text {
    background-color: #f8f9fa;
    padding: 1rem;
    border-radius: 5px;
    margin: 1rem 0;
    font-size: 0.9rem;
    color: #495057;
}

.overview-card {
    background-color: #ffffff;
    border-radius: 10px;
    padding: 2rem;
    margin: 1rem 0;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
}

.card-title {
    color: #2c3e50;
    border-bottom: 2px solid #e9ecef;
    padding-bottom: 0.5rem;
    margin-bottom: 1rem;
}

.info-card {
    background-color: #ffffff;
    border-radius: 10px;
    padding: 2rem;
    margin: 1rem 0;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
}

.info-text {
    line-height: 1.6;
    color: #495057;
}

.feature-list {
    list-style-type: none;
    padding-left: 0;
    margin-top: 1rem;
}

.feature-list li {
    margin-bottom: 0.5rem;
    padding-left: 1.5rem;
    position: relative;
}

.feature-list li:before {
    content: "•";
    color: #2c3e50;
    position: absolute;
    left: 0;
}

.required-field {
    color: #dc3545;
    font-size: 0.8rem;
    margin-top: 0.25rem;
}

.centered-container {
    color: #2c3e50;
    border-bottom: 2px solid #e9ecef;
    padding-bottom: 0.5rem;
    margin-bottom: 1rem;
    max-width: 66%;
    margin: 0 auto;
}

.upload-instruction {
    margin-bottom: 0;
    line-height: 1.2;
}

.required-field {
    color: #dc3545;
    margin-top: 0.2rem;
}

:root {
    --primary-color: #2C3E50;
    --secondary-color: #3498DB;
    --accent-color: #E74C3C;
    --background-color: #ECF0F1;  /* Change this line */
    --text-color: #2C3E50;
    --card-background: #FFFFFF;
}

.scrollable-text {
    height: 300px;  /* Adjust this value to control the height of the scrollable area */
    overflow-y: auto;
    padding: 15px;
    border: 1px solid #e0e0e0;
    border-radius: 5px;
    background-color: rgba(255, 255, 255, 0.8);
}

.transparent-card {
        background-color: rgba(255, 255, 255, 0.2);  /* White with 20% opacity */
        backdrop-filter: blur(10px);  /* Optional: adds a blur effect */
        border-radius: 10px;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        margin-bottom: 2rem;
        border: none;
        transition: transform 0.2s;
    }

    /* Optional hover effect */
    .transparent-card:hover {
        background-color: rgba(255, 255, 255, 0.3);  /* Slightly more opaque on hover */
    }
.navbar {
        position: sticky;
        top: 0;
        z-index: 1000;  /* Ensures navbar stays on top of other content */
        background-color: var(--surface-color);  /* Use your navbar color */
    }
"""